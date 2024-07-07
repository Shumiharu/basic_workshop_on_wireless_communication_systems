///////// OFDM-QPSKのAWGN通信路におけるBER特性  //////////////////

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1						//波数
#define s 8						//オーバーサンプル比
#define K 1					//サブキャリア数
#define N (s*K)					//全サンプル数
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1シンボル時間
#define fc a/Ts 				//搬送波周波数fc
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

//BER特性用↓↓↓↓
#define bit 2.0 				//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1   				
#define loop pow(10,4)			//試行回数

////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void QPSK(double complex base[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void SP(double complex Serial[],double complex Parallel[]);
void PS(double complex Parallel[],double complex Serial[]);


int main(void)
{	
	int i,j,m; 
	
///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 

////////////////////// 時間と周波数の設定　　//////////////////////
	double t[N]; 						//時間t    0～K*Tsまでの時間
	double f[N];						//周波数f  -Ts/s～Ts/sまでの周波数
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}
	
	double complex sF1_QPSK[N];		//QPSK信号
	double complex sF2_QPSK[N];		//先頭サンプルのみQPSK信号
	double complex sF[N];			//S/P変換後の信号
	double complex st[N];			//IFFT後の信号	
	double complex awgn[N];			//ガウス雑音
	double complex rt[N];			//
	double complex rF[N];			//
	double complex ro[N];			//
	
	QPSK(sF1_QPSK);					//QPSK信号生成
	inter(sF1_QPSK,sF2_QPSK);		//K本のサブキャリア=先頭サンプルのみ抜き出す
	SP(sF2_QPSK,sF);
	IFFT(sF,st);
	
	// putchar('\n');
	// for(i=0;i<N;i++)
		// printf("%f\t%f\t%f\n",t[i],creal(st[i]),cimag(st[i]));
	
///////////////////////////////////BER特性算出時に用いる変数　　/////////////////////////////////////////////////
	double sigma;				//ガウス雑音の分散			
	int EbN0;					//Eb/N0
	double ci,cq;				//BERカウント用
	double ber;  				//誤り率
	
	printf("s=%d\n",s);
	printf(" ρ[dB]         BER \n");
	for(EbN0=0; EbN0<=12; EbN0++){
		ci=0.0; cq=0.0; 
		sigma=(double)BnTs/(2*s*bit)*pow(10.0,(-0.1*EbN0));		//ガウス雑音の分散算出
		for(m=0; m<loop; m++){	
			///////　ガウス雑音付加　　////////////////////////////////////////////////////////////////////////
			boxmuller(N,awgn,0.0,sigma);		//　AWGNの生成　
			for(i=0;i<N;i++)
				rt[i]=st[i]+awgn[i];		//rt(t)=st(t)+AWGN(t)　
			
			///////  FFT r(t)→R(f) ///////////////////////////////////////////////////////////////////
			FFT(rt,rF);
			
			//////  P/S変換　　///////////////////////////////////////////////////////////////////////////  
			PS(rF,ro);
			
			///////// 各信号の時間波形とスペクトルの確認　//////////////////////////////////////////////////////
			/* putchar('\n');
			printf("Time\ts_I(t)\tr_I(t)\n");
			for(i=0;i<N;i++)
				printf("%f\t%f\t%f\n",t[i],creal(st[i]),creal(rF[i]));
			
			putchar('\n');
			printf("Freq\t|s(f)|\t|r(f)|\n");
			for(i=0;i<N;i++)
				printf("%f\t%e\t%e\n",f[i],cabs(sF[i]),cabs(rF[i]));
			
			putchar('\n');
			for(i=0;i<N;i++)
				printf("%f\t%f\t%f\t%e\n",t[i],creal(ro[i]),f[i],cabs(ro[i]));
			*/
			
			
			//// 誤り判定(先頭点のみ)  ////////////////////////////////////////////////////////////////////////
			for(j=0; j<K; j++){
				if(creal(sF2_QPSK[j*s])>=0 && creal(ro[j*s])<0) ci++;
				if(creal(sF2_QPSK[j*s])<0 && creal(ro[j*s])>=0) ci++;
				if(cimag(sF2_QPSK[j*s])>=0 && cimag(ro[j*s])<0) cq++;	
				if(cimag(sF2_QPSK[j*s])<0 && cimag(ro[j*s])>=0) cq++;
			}	
		}
	ber=(ci+cq)/(bit*(double)loop*K);
	printf("%d\t%e\n",EbN0,ber);
	}
	

	
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
void QPSK(double complex base[]){
		int i,j;
		double tmp_I;	//同相成分Iを一時的に格納
		double tmp_Q;	//直交成分Qを一時的に格納
		
		for(i=0; i<K; i++){
			tmp_I=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
			tmp_Q=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
			
			for(j=0; j<s; j++)
				base[s*i+j]=tmp_I+tmp_Q*I;				//複素ベースバンド信号u(t)=I(t)+jQ(t)
		}			
}

////////////////　先頭サンプルのみ抜き出し関数　/////////////////////////
void inter(double complex s1[],double complex s2[]){
	int i,j;
	
	for(i=0;i<N;i+=s){		//シンボル数だけloop
			s2[i]=s1[i];		
			for(j=1;j<s;j++)	//1サンプル以降を0にする
				s2[i+j]=0.0+0.0*I;
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
/////////////////////////Serial to Pallel/////////////////////////
void SP(double complex Serial[],double complex Parallel[]){
	int i,j;
	
	j=0;
	for(i=0; i<N; i++){
		if(i>=(N-K)/2 && i<(N+K)/2){
			Parallel[i]=Serial[j*s];
			j++;
		}
		else{
			Parallel[i]=0.0+0.0*I;
		}
	}	
}

////////////////////////Parallel to Serial///////////////////////
void PS(double complex Parallel[],double complex Serial[]){
	int i,j;
	
	j=0;
	for(i=0; i<N; i++){
		if(i>=(N-K)/2 && i<(N+K)/2){
			Serial[j*s]=Parallel[i];
			j++;
		}
	}	
}


	










