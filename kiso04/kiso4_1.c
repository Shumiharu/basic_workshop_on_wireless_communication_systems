////////AWGN通信路のルートロールオフ伝送系におけるQPSK伝送のBER特性/////////////////
//・帯域制限
#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム

#define Ts 0.001				//1シンボル時間，サンプリング周期
#define fs 1/Ts					//サンプリング周波数
#define a 1						//波数
#define s 1					//オーバーサンプル比
#define K 256	//シンボル数
#define f0 (a/Ts)				//変調周波数
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define N (s*K)			//全サンプル数
#define delta_f (1/(Ts*K))	//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長
#define loop 10000			//試行回数
//BER特性用↓↓↓↓
#define c 2 					//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1    				//1シンボル当たりb点を用いるとき，BnTs=b

/////////////////////////　　　　FFT　　　　//////////////////////////////
void FFT(double complex Time[] , double complex Freq[]){

	int i,m,k,j;
	double arg;
	double complex X[N],Y[N];
	
	for(i=0;i<N; i++){
		Freq[i]=0;
		X[i]=Time[i];
	}
	for(m=2; m<=N; m*=2){		//バタフライの段数
		j=m/2;		//何個とびかを表す
		
		for(i=0; i<N/m; i++){		//
			for(k=0; k<j; k++){			//どこの数字と計算するか，何個とびか
				arg=2*M_PI*(k-N/2)/m;
				Y[i*m+k]=X[i*j+k]+cexp(-arg*I)*X[i*j+(N/2)+k];
				Y[i*m+j+k]=X[i*j+k]-cexp(-arg*I)*X[i*j+(N/2)+k];
			}
		}
		for(i=0; i<N; i++){
			X[i]=Y[i];
		}
		for(i=0; i<N; i++)
			Y[i]=0.0+0.0*I;
	}
	
	for(i=0;i<N;i++)
		Freq[i]=X[i]/sqrt((double)N);
}

void IFFT(double complex Freq[],double complex Time[]){
	int i,j,k,m;
	double arg;
	double complex X[N],Y[N];
	
	for(i=0;i<N/2;i++){
		X[i]=Freq[i+N/2];
		X[i+N/2]=Freq[i];
		Time[i]=0;
		Time[i+N/2];
	}
	
	
	for(m=2;m<=N;m*=2){
		j=m/2;
		for(i=0;i<N/m;i++){
			for(k=0;k<j;k++){
				arg=2*M_PI*(k)/m;
				Y[i*m+k]=X[i*j+k]+cexp(arg*I)*X[i*j+(N/2)+k];//バタフライ上側
				Y[i*m+j+k]=X[i*j+k]-cexp(arg*I)*X[i*j+(N/2)+k];//バタフライ下側
			}
		}
		for(i=0;i<N;i++){
			X[i]=Y[i];
		}
		for(i=0;i<N;i++){
			Y[i]=0.0+0.0*I;
		}
	}
	for(i=0;i<N;i++){
		Time[i]=X[i]/sqrt(N);
	}
}

//離散フーリエ逆変換(Inverse Disecrete Fourier Transform)を行う関数//	
void idft(double complex spectrum[] , double complex signal[]){
	for(int n=0; n<=N-1; n++){
		signal[n]=0;
	}
	
	for(int n=0; n<=N-1; n++){
		for(int k=0; k<=N-1; k++){
			signal[n]+=(1/sqrt((double)N))*spectrum[k]*cexp(I*2*M_PI*(k-N/2)*(n-N/2)/N);
		//　　本当はkやnに-N/2～を入れたいが，C言語の配列は0からなのでk-N/2やn-N/2にする　　//
		}
	}	
}


/////////////   複素ベースバンド信号(と搬送帯域信号)の生成　　　//////////////////////
///　　搬送帯域信号も欲しいなら☆の部分のコメント化を解除する   ////
void baseband(const double time[] , double complex base[]){
		double tmp_I;	//同相成分Iを一時的に格納
		double tmp_Q;	//直交成分Qを一時的に格納
		
		for(int i=0; i<K; i++){
			tmp_I=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
	
			tmp_Q=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
			for(int j=0; j<s; j++){
				base[s*i+j]=tmp_I+tmp_Q*I;				//複素ベースバンド信号u(t)=I(t)+jQ(t)
				// printf("%f  ", time[s*i+j]);
				// printf("%f  ", creal(base[s*i+j]));
				// printf("%f  ", cimag(base[s*i+j]));
			}
		}	
}

////////　　　　Box-Muller法によるガウス雑音生成　　　　　　///////
void boxmuller(double complex noise[], double mean , double variance){
	double u1,u2;
	double nreal[N],nimag[N];
	
	for(int i=0; i<N; i++){
	////////////////　　　　[0,1]の一様乱数生成　　　　　/////////////

		u1=((double)rand()+1.0)/((double)RAND_MAX+1.0);		//rand()=0～32767
		u2=((double)rand()+1.0)/((double)RAND_MAX+1.0);
		
		///   標準偏差σを1/√2にすれば各軸の分散σ^2=1/2になる．
		nreal[i]=sqrt(variance)*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2)+mean;		//分散を1/√2にして合成後の電力を1Wにする．
		nimag[i]=sqrt(variance)*sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2)+mean;
		noise[i]=nreal[i]+nimag[i]*I;
	}
}


int main(void)
{	
	int i,j;  
	
	
	
	////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 
	
	/////////////// 時間と周波数の設定　　/////////////////////
	double t[N]; 						//時間t    0～K*Tsまでの時間
	double f[N];						//周波数f  -Ts/s～Ts/sまでの周波数
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}
	
	//////////////  ロールオフフィルタのインパルス応答(図4-1)  //////////
	double complex h1[N];	//ルートロールオフフィルタのインパルス応答
	double alpha1=0.5;		//ロールオフ率（Roll-off Ratio）
	double time1[N],arg1;
	
	for(i=0;i<N;i++){
		time1[i]=T_sample*(i-N/2);
		arg1=M_PI*fs*(time1[i]);
			
		if(time1[i]==0.0)
			h1[i]=sqrt((double)N);		// 0/0となる部分の極限は0だから
		else if(time1[i]==1/(2*alpha1*fs) || time1[i]==-1/(2*alpha1*fs))	
			h1[i]=sqrt((double)N)*alpha1/2*sin(M_PI/(2*alpha1));		// 0/0となる部分の極限はπ/4だから
		else
			h1[i]=sqrt((double)N)*sin(arg1)/arg1*cos(alpha1*arg1)/(1-pow(2*alpha1*fs*time1[i],2));
	}
	
	
	
	////////////　　ロールオフフィルタの周波数応答と結果表示　　////////////////
	double complex fft_h1[N];			//ルートロールオフフィルタのh1(t)の周波数成分
	
	FFT(h1,fft_h1);
	// putchar('\n');
	// printf("　時間　　　　　インパルス応答　　　　   周波数　　　　　 　周波数応答 \n");
	// for(i=0;i<N;i++)
		// printf(" %f          %f           %f            %f \n",time1[i],creal(h1[i]),f[i],cabs(fft_h1[i]));
		
	//////////////　　フィルタリング前後の時間波形と周波数スペクトル　///////////////
	double complex s1[N];				//ベースバンド信号
	double complex fft_s1[N];			//ベースバンド信号のスペクトル
	double complex s2[N];				//フィルタ通過後の時間信号
	double complex fft_s2[N];			//フィルタ通過後の周波数スペクトル
	double complex s11[N];				//ベースバンド信号の1サンプル目，Tsごとの離散信号にするため
	double complex fft_s11[N];			//Ts毎の離散信号のスペクトル
	
	baseband(t,s1);		//ベースバンド信号生成
	
	for(i=0;i<N;i+=s){		//シンボル数だけloop
		s11[i]=s1[i];		
		for(j=1;j<s;j++)	//1サンプル以降を0にする
			s11[i+j]=0.0+0.0*I;
	}
	// for(i=0;i<N;i++){
		// printf("%f  %f   %f    %f  \n",creal(s1[i]),cimag(s1[i]),creal(s11[i]),cimag(s11[i]));
	// }確認用
	
	FFT(s1,fft_s1);		//ベースバンド信号のスペクトル
	FFT(s11,fft_s11);	//Tsごとのベースバンド信号のスペクトル
	
	for(i=0;i<N;i++){
		fft_s2[i]=fft_s11[i]*fft_h1[i];		//S2(f)=s11(f)*H1(f)
	}
	IFFT(fft_s2,s2);		//s2(t)=IDFT{S2(f)}離散フーリエ逆変換で時間信号に戻す	
	
	////////////////　フィルタリング前後の時間，周波数波形の表示 ///////////////
/* 	putchar('\n');
	printf(" 　時間　　　　 　フィルタ前s1_I　　   フィルタ後s2_I      フィルタ前s1_Q      フィルタ後s2_Q\n");
	for(i=0;i<N;i++)
		printf(" %f          %f            %f           %f            %f\n",t[i],creal(s1[i]),creal(s2[i]),cimag(s1[i]),cimag(s2[i]));
		 */
	
/* 	putchar('\n');
	printf(" 　周波数　　　　 　 フィルタ前|S1(f)|　　　　   フィルタ後|S2(f)| \n");	
	for(i=0;i<N;i++)
		printf(" %f          %e            %e\n",f[i],cabs(fft_s1[i]),cabs(fft_s2[i]));
	 */
	/////////////////　ルートロールオフ伝送系  ////////////////////////
	double complex gt[N],fft_gt[N];			//送信フィルタの時間波形gt(t),周波数成分Gt(f)
	double complex gr[N],fft_gr[N];			//送信フィルタの時間波形gr(t),周波数成分Gr(f)
	double complex s3[N],fft_s3[N];			//送信フィルタに入力するQPSK信号の時間波形s3(t)，周波数スペクトルS3(f)
	double complex s4[N],fft_s4[N];			//送信フィルタ通過後の時間信号s4(t)，周波数スペクトルS4(f)
	double complex awgn[N],fft_awgn[N];		//AWGNの時間波形とスペクトル
	double complex r0[N],fft_r0[N];			//AWGNが加わった時間信号r0(t),周波数スペクトルR0(f)
	double complex r[N],fft_r[N];			//受信フィルタ通過後の時間信号r(t),周波数スペクトルR(f)
	double sigma;							//分散			
	int EbN0;								//Eb/N0
	double ci,cq;							//BERカウント用
	double ber;  							//誤り率
	double attempt;
	
	//////////////　ルート配分　///////////////
	for(i=0;i<N;i++){
		fft_gt[i]=sqrt(cabs(fft_h1[i]))+I*0;	//sqrt(double)なので実部，虚部分けてやる
		fft_gr[i]=fft_gt[i];
	}
	
	//////////////  ベースバンド信号生成　////////////////
		baseband(t,s1);		//ベースバンド信号生成
			
		for(i=0;i<N;i+=s){		//シンボル数だけloop
			s11[i]=s1[i];		
			for(j=1;j<s;j++)	//1サンプル以降を0にする
				s11[i+j]=0.0+0.0*I;
		}
	
		////////////　QPSK複素ベースバンド信号s0(t)生成　周波数スペクトルS0(f)　/////////////////
		for(i=0;i<N;i++)
			s3[i]=s11[i];
			
		FFT(s3,fft_s3);		//s3(t)→S3(f)
			

		////////////  送信フィルタ通過　　///////////////
		for(i=0;i<N;i++)		
			fft_s4[i]=fft_gt[i]*fft_s3[i];		//S4(f)=Gt(f)*S3(f)
			
		IFFT(fft_s4,s4);		//S4(f)→s4(t)	
	
	
	printf("α=%.1f  s=%d   \n",alpha1,s);
	printf(" ρ[dB]         BER \n");
	for(EbN0=0; EbN0<=9; EbN0++){
		ci=0.0; cq=0.0;
		attempt=0.0;
		sigma=(double)BnTs/(2*s*c)*pow(10.0,(-0.1*EbN0));
		for(int m=0; m<loop; m++){
			////////////　AWGNの生成　　//////////////////	
			boxmuller(awgn,0.0,sigma);		//AWGN生成
			
			///////////  s4(t)+AWGN(t)　/////////////////
			for(i=0;i<N;i++)
				r0[i]=s4[i]+awgn[i];		//雑音が加わる r0(t)=s4(t)+awgn(t)
			
			///////////  受信フィルタ通過  /////////////
			FFT(r0,fft_r0);			//r0(t)→R0(f)
			
			for(i=0;i<N;i++)
				fft_r[i]=fft_r0[i]*fft_gr[i];	//受信フィルタR(f)=R0(f)*Gr(f)
			
			IFFT(fft_r,r);		//受信フィルタ通過後の時間波形r(t)
			
			///////////////　　信号の時間波形，スペクトルを表示　　///////////////
			////////////////　　各時間信号波形の表示（I軸）　　///////////////////////
			// putchar('\n');
			// printf("　時間　　　       s0_I             s_I               r0_I                r_I \n");
			// for(i=0;i<N;i++)
				// printf(" %f       %f         %f          %f           %f\n",t[i],creal(s1[i]),creal(s4[i]),creal(r0[i]),creal(r[i]));
			
			////////////////　　各時間信号波形の表示（Q軸)  ////////////////////////
			// putchar('\n');
			// printf("　時間　　　       s0_Q             s_Q               r0_Q                r_Q \n");
			// for(i=0;i<N;i++)
				// printf(" %f       %f         %f          %f           %f\n",t[i],cimag(s1[i]),cimag(s4[i]),cimag(r0[i]),cimag(r[i]));
			
			/////////////////　　各スペクトルの表示　　　///////////////////////
			//                     送信フィルタ通過後　　　　　　　　　AWGN加算後　　　　　　　　　　受信フィルタ通過後
			// putchar('\n');
			// printf("　周波数　　　           |S(f)|                 |R0(f)|                  |R(f)|  \n");
			// for(i=0;i<N;i++)
				// printf(" %f          %e           %e            %e \n",f[i],cabs(fft_s4[i]),cabs(fft_r0[i]),cabs(fft_r[i]));
			
			for(int j=0; j<K; j++){
				////////////  BER  ////////////
				if(creal(s1[j*s])>=0 && creal(r[j*s])<0) ci++;
				if(creal(s1[j*s])<0 && creal(r[j*s])>=0) ci++;
				if(cimag(s1[j*s])>=0 && cimag(r[j*s])<0) cq++;	
				if(cimag(s1[j*s])<0 && cimag(r[j*s])>=0) cq++; 
				
	/* 			////////////  SER   ///////////
				if(creal(s1[j*s])>=0 && cimag(s1[j*s])>=0){	//信号点が第一象限の時
					if(creal(r[j*s])<0 || cimag(r[j*s])<0) cs++;
				}
				if(creal(s1[j*s])<0 && cimag(s1[j*s])>=0){	//信号点が第二象限の時
					if(creal(r[j*s])>=0 || cimag(r[j*s])<0) cs++;
				}
				if(creal(s1[j*s])<0 && cimag(s1[j*s])<0){		//信号点が第三象限の時
					if(creal(r[j*s])>=0 || cimag(r[j*s])>=0) cs++;
				}
				if(creal(s1[j*s])>=0 && cimag(s1[j*s])<0){	//信号点が第四象限の時
					if(creal(r[j*s])<0 || cimag(r[j*s])>=0) cs++;
				} */
			}	
		attempt++;
		}	
	ber=(ci+cq)/(2.0*(double)loop*K);
	printf(" %d          %e \n",EbN0,ber);
	}	
	
	

	return 0;
}





