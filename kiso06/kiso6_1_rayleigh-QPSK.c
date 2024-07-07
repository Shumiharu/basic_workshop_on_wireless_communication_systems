////////// 動的または準静的レイリーフェージング環境下でのQPSKのBER特性，時間領域等化 /////////////
///////ルートロールオフ伝送系///////////

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1						//波数
#define s 4						//オーバーサンプル比
#define K 1024					//1フレーム長
#define N (s*K)					//全サンプル数
#define frame 100				//フレーム数
#define sigmaR 1.0				//レイリー分布の分散
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1シンボル時間
#define fs 1/Ts					//サンプリング周波数
#define fc a/Ts 				//搬送波周波数fc
#define fdTs 0.01				//fdTs
#define fd (fdTs/Ts)			//ドップラーシフト
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

#define loop pow(10,7)					//試行回数

//BER特性用↓↓↓↓
#define c 2.0 					//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1   				

////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void baseband(const double time[] , double complex base[]);
void boxmuller(int num,double complex noise[],double mean,double variance);
void Rayleigh(double complex R[]);		//動的レイリーフェジングの生成関数
void qs_Rayleigh(double complex qsRfading[]);		//準静的フェージングの生成関数

int main(void)
{	
	int i,j,m; 
	
///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 
	
	//ファイルへの書き込み
	FILE *fp;
	char *fname="result_kiso6_1.csv";
	fp=fopen(fname, "w");
	if(fp==NULL){
	printf("%sファイルが開けません\n", fname);
		return (-1);
	}
	else{
	
//////////////////////// 時間波形，スペクトル表示用の時間と周波数の設定　　/////////////////////
	double *t; 						//時間t
	double *f;						//周波数f

	t=(double *)malloc(sizeof(double)*N);		//メモリ拡張
	f=(double *)malloc(sizeof(double)*N);

	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}	
	
////////////////////////　　フェージングの作成に用いる変数　　/////////////////////////
	
	

///////////////////////////////　ルートロールオフ伝送系構築に用いる変数  ////////////////////////
	double complex *s3;					//送信フィルタに入力するQPSK信号の時間波形s3(t)
	double complex *s31,*fft_s31;		//QPSK信号の先頭サンプルのみの時間波形とそのスペクトル
	double complex *s4,*fft_s4;			//送信フィルタ通過後の時間信号s4(t)，周波数スペクトルS4(f)
	double complex *s5;					//動的レイリーフェージングを受けた時間信号s5(t)
	double complex *awgn;				//AWGNの時間波形
	double complex *Rfading;			//動的または準静的レイリーフェージング
	double complex *r0,*fft_r0;			//AWGNが加わった時間信号r0(t),周波数スペクトルR0(f)
	double complex *r,*fft_r;			//受信フィルタ通過後の時間信号r(t),周波数スペクトルR(f)
	double complex *eq_r;				//等化後の受信信号r(t)
	
	h1=(double complex *)malloc(sizeof(double complex)*N);
	time1=(double *)malloc(sizeof(double)*N);
	fft_h1=(double complex *)malloc(sizeof(double complex)*N);
	gt=(double complex *)malloc(sizeof(double complex)*N);
	fft_gt=(double complex *)malloc(sizeof(double complex)*N);
	gr=(double complex *)malloc(sizeof(double complex)*N);
	fft_gr=(double complex *)malloc(sizeof(double complex)*N);
	s3=(double complex *)malloc(sizeof(double complex)*N);
	s31=(double complex *)malloc(sizeof(double complex)*N);
	fft_s31=(double complex *)malloc(sizeof(double complex)*N);
	s4=(double complex *)malloc(sizeof(double complex)*N);
	fft_s4=(double complex *)malloc(sizeof(double complex)*N);
	s5=(double complex *)malloc(sizeof(double complex)*N);
	awgn=(double complex *)malloc(sizeof(double complex)*N);
	Rfading=(double complex *)malloc(sizeof(double complex)*N);
	r0=(double complex *)malloc(sizeof(double complex)*N);
	fft_r0=(double complex *)malloc(sizeof(double complex)*N);
	r=(double complex *)malloc(sizeof(double complex)*N);
	fft_r=(double complex *)malloc(sizeof(double complex)*N);
	eq_r=(double complex *)malloc(sizeof(double complex)*N);
	
/////////////////////////////////////　ロールオフフィルタのインパルス応答算出　　//////////////////////////////
	double complex *h1;						//ルートロールオフフィルタのインパルス応答
	double alpha1=0.5;						//ロールオフ率（Roll-off Ratio）
	double *time1,arg1;
	double complex *fft_h1;				//ルートロールオフフィルタのh1(t)の周波数成分
		
	double complex *fft_gt;			//送信フィルタの時間波形gt(t),周波数成分Gt(f)
	double complex *fft_gr;			//送信フィルタの時間波形gr(t),周波数成分Gr(f)
	
	h1=(double complex *)malloc(sizeof(double complex)*N);
	time1=(double *)malloc(sizeof(double)*N);
	fft_h1=(double complex *)malloc(sizeof(double complex)*N);
	fft_gt=(double complex *)malloc(sizeof(double complex)*N);
	fft_gr=(double complex *)malloc(sizeof(double complex)*N);
	
	
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
//////////////////////////////////　　ロールオフフィルタの周波数応答　　//////////////////////////////////////////
	FFT(h1,fft_h1);
	
///////////////////////////////////////////　ルート配分　////////////////////////////////////////////////
	for(i=0;i<N;i++){
		fft_gt[i]=sqrt(cabs(fft_h1[i]))+I*0;	//sqrt(double)なので実部，虚部分けてやる
		fft_gr[i]=fft_gt[i];	
	}
	
////////////////////////////////////////////  ベースバンド信号生成　////////////////////////////////////////////
		baseband(t,s3);		//ベースバンド信号生成
			
		for(i=0;i<N;i+=s){		//シンボル数だけloop
			s31[i]=s3[i];		
			for(j=1;j<s;j++)	//1サンプル以降を0にする
				s31[i+j]=0.0+0.0*I;
		}
		FFT(s31,fft_s31);		//s31(t)→S31(f)
		
//////////////////////////////////////////  送信フィルタ通過s31(t)→s4(t)　　////////////////////////////////////////////////
		for(i=0;i<N;i++)		
			fft_s4[i]=fft_gt[i]*fft_s31[i];		//S4(f)=Gt(f)*S31(f)
			
		IFFT(fft_s4,s4);		//S4(f)→s4(t)	

////////////////////////////////////BER特性算出時に用いる変数　　/////////////////////////////////////////////////
	double sigma;				//ガウス雑音の分散			
	int EbN0;					//Eb/N0
	double ci,cq;				//BERカウント用
	double ber;  				//誤り率
	
	printf("動的または準静的フェージング環境下におけるQPSK伝送のBER特性\n");
	printf("s=%d  fdTs=%f \n",s,fdTs);
	printf(" ρ[dB]         BER \n");
	
	for(EbN0=0; EbN0<=0; EbN0+=5){
		ci=0.0; cq=0.0; 
		sigma=(double)BnTs/(2*s*c)*pow(10.0,(-0.1*EbN0));		//ガウス雑音の分散算出
		// do{
		for(m=0; m<1; m++){
			
			
			Rayleigh(Rfading);					//　動的レイリーフェージングの生成
			// qs_Rayleigh(Rfading);			//準静的レイリーフェージングの生成
			
			boxmuller(N,awgn,0.0,sigma);		//　AWGNの生成　　			
			
			///// 動的または準静的レイリーフェージング付加　////////////////////
			for(i=0;i<N;i++)
				s5[i]=s4[i]*Rfading[i];		// s5(t)=s4(t)*Rfading(t)  動的レイリーフェージングが加わる　
			
			
			////////  動的または準静的レイリーフェージングによる時間波形とスペクトルの変化  /////////
		/* 	double complex fft_s4[N],fft_s5[N];
			FFT(s4,fft_s4);
			FFT(s5,fft_s5);
			
			// fprintf(fp,"TIme\ts4(t)\ts5(t)\n");
			fprintf(fp,"Freq\t|S4(f)|\t|S5(f)|\n");
			for(i=0;i<N;i++){
				// fprintf(fp,"%f\t%f\t%f\n",t[i],s4[i],s5[i]);
				fprintf(fp,"%f\t%e\t%e\n",f[i],cabs(fft_s4[i]),cabs(fft_s5[i]));
			 } */
			
			///////　ガウス雑音付加　　//////////////////////
			for(i=0;i<N;i++){
				r0[i]=s5[i]+awgn[i];		//r0(t)=s5(t)+AWGN(t)　
			}			
			////////　受信フィルタ通過  /////////////////////
			FFT(r0,fft_r0);		//r0(t)→R0(f)
			
			for(i=0;i<N;i++)
				fft_r[i]=fft_r0[i]*fft_gr[i];	//受信フィルタ通過R(f)=R0(f)*Gr(f)
			
			IFFT(fft_r,r);		//受信フィルタ通過後の時間波形r(t)	
	
			/////////////////////　 動的または準静的レイリーフェージング等化  /////////////////////////////
			for(i=0;i<N;i++)
				r[i]=r[i]/Rfading[i];
		
			///////////////////////////  動的または準静的レイリーフェージング誤り判定(先頭点のみ)  /////////////////////////
			for(j=0; j<K; j++){
				if(creal(s3[j*s])>=0 && creal(r[j*s])<0) ci++;
				if(creal(s3[j*s])<0 && creal(r[j*s])>=0) ci++;
				if(cimag(s3[j*s])>=0 && cimag(r[j*s])<0) cq++;	
				if(cimag(s3[j*s])<0 && cimag(r[j*s])>=0) cq++;
			}	
		}
	ber=(ci+cq)/(c*(double)loop*K);
	printf("%d\t%e\n",EbN0,ber);
	}

	
	
	
	}
	
	fclose(fp);
	printf("%s 書き込み完了\n", fname);
	
	
	return 0;
}

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

////////////  高速フーリエ逆変換IFFT(Inverse FFT)　　/////////////
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

/////////////   複素ベースバンド信号(と搬送帯域信号)の生成　　　//////////////////////
void baseband(const double time[] , double complex base[]){
		int i,j;
		double tmp_I;	//同相成分Iを一時的に格納
		double tmp_Q;	//直交成分Qを一時的に格納
		
		for(i=0; i<K; i++){
			tmp_I=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
	
			tmp_Q=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
			for(j=0; j<s; j++){
				base[s*i+j]=tmp_I+tmp_Q*I;				//複素ベースバンド信号u(t)=I(t)+jQ(t)
				// printf("%f  ", time[s*i+j]);
				// printf("%f  ", creal(base[s*i+j]));
				// printf("%f  ", cimag(base[s*i+j]));
			}
		}	
}



//　　　　Box-Muller法によるnum個からなる配列=ガウス雑音生成　　　　　　//
void boxmuller(int num,double complex noise[], double mean , double variance){
	int i,j;
	double u1,u2;
	double nreal[num],nimag[num];
	
	for(i=0; i<num; i++){
	////////////////　　　　[0,1]の一様乱数生成　　　　　/////////////

		u1=((double)rand()+1.0)/((double)RAND_MAX+1.0);		//rand()=0～32767
		u2=((double)rand()+1.0)/((double)RAND_MAX+1.0);
		
		///   標準偏差σを1/√2にすれば各軸の分散σ^2=1/2になる．
		nreal[i]=sqrt(variance)*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2)+mean;		//分散を1/√2にして合成後の電力を1Wにする．
		nimag[i]=sqrt(variance)*sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2)+mean;
		noise[i]=nreal[i]+nimag[i]*I;
	}
}
///////////////////////////// 動的レイリーフェージングの生成　　///////////////////////////////////
void Rayleigh(double complex R[]){
	int i,j;
	double *freq;						//周波数
	double complex *noise;				//ガウス雑音
	double complex *fft_noise;			//ガウス雑音のスペクトル
	// double complex *Rfading;			//レイリーフェージングの時間波形
	double complex *fft_R;		//レイリーフェージングのスペクトル
	double complex *F_tmp;				//(レイリーフェージングのスペクトル)×(ガウス雑音のスペクトル)
	double *Pf_Rfading;					//レイリーフェージングの電力スペクトル密度
	double PP_ray=0.0;
	////////////// マロック関数による領域拡張　//////////////////
	freq=(double *)malloc(sizeof(double)*N);
	noise=(double complex *)malloc(sizeof(double complex)*N);	
	fft_noise=(double complex *)malloc(sizeof(double complex)*N);
	// Rfading=(double complex *)malloc(sizeof(double complex)*N);
	fft_R=(double complex *)malloc(sizeof(double complex)*N);
	F_tmp=(double complex *)malloc(sizeof(double complex)*N);
	Pf_Rfading=(double *)malloc(sizeof(double)*N);
	
///////////////////////////// 周波数の計算　　/////////////////////////////
	for(i=0; i<N; i++)
		freq[i]=(i-N/2)*delta_f;
	
//////////////////////////////　電力スペクトル密度算出　/////////////////////////////////
	for(i=0;i<N;i++){
		if(fd>fabs(freq[i]))	//sqrt()の中身が負にならないように
			Pf_Rfading[i]=sigmaR/M_PI/(sqrt(fd*fd-freq[i]*freq[i]));
		else 
			Pf_Rfading[i]=0.0;
		
		PP_ray+=Pf_Rfading[i];
	}
///////////////////////////////　ガウス雑音の周波数成分 //////////////////////////	
	boxmuller(N,noise,0.0,sigmaN);
	FFT(noise,fft_noise);
	
	for(i=0;i<N;i++){
		fft_R[i]=sqrt(Pf_Rfading[i]*N/PP_ray);		//レイリーフェージングのスペクトル？？？？
		F_tmp[i]=fft_noise[i]*fft_R[i];		//レイリーフェージングのスペクトル×AWGNのスペクトル
	}
	
///////　(レイリーフェージングのスペクトル)×(ガウス雑音のスペクトル)を時間領域に戻す=フェージング完成  ///////
	IFFT(F_tmp,R);
}

void qs_Rayleigh(double complex qsRfading[]){
	int i,j;
	double complex noise[1];		//ガウス雑音

		////////////////////　　ガウス雑音の生成　　////////////////////////
		boxmuller(1,noise,0.0,sigmaN);		//ガウス雑音の生成
			
		////////////　　準静的レイリーフェージングの作成＝ガウス雑音を1フレームずつコピーする  ////////////////
		for(i=0;i<N;i++){	//1フレーム単位でコピー
			qsRfading[i]=noise[0];		//1フレームにN=s*K点すべてに同じawgnを入れる．
		}
}

/* 
void qs_Rayleigh(double complex qsRfading[]){
	int i,m;
	double complex *noise;		//ガウス雑音
	double *t2;					//準静的レイリーフェージング用の時間
		////////////// マロック関数による領域拡張　//////////////////
		noise=(double complex *)malloc(sizeof(double complex)*frame);		//配列の要素数をframeにしてもできる　m=0～frame-1までの計frame個だから			
		
		////////////////////　　ガウス雑音の生成　　////////////////////////
		boxmuller(frame,noise,0.0,sigmaN);		//ガウス雑音の生成
			
		for(m=0;m<frame;m++){		//frame数繰り返す	
			////////////　　準静的レイリーフェージングの作成＝ガウス雑音を1フレームずつコピーする  ////////////////
			for(i=0;i<K;i++){	//1フレーム単位でコピー
				qsRfading[i+m*K]=noise[m];		//1フレーム，つまりK16点分同じガウス雑音をコピーする，awgnはframe分配列があればOK
				// t2[i+m*K]=Ts*(i+m*K);		//準静的フェージング1フレーム(K=16点)でK*Ts秒，つまり1点でTs秒割り振る
			}
		}
}

 */



