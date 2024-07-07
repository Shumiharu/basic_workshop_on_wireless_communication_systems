//////////動的レイリーフェージングの振幅，位相波形及び電力の確認

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1						//波数
#define s 8						//オーバーサンプル比
#define K 1024					//1フレーム長
#define N (s*K)					//全サンプル数

#define Ts 0.001				//1シンボル時間
#define fc a/Ts 				//搬送波周波数fc
#define fdTs 0.01				//fdTs
#define fd (fdTs/Ts)			//ドップラーシフト
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

#define loop 10					//試行回数
//ガウス雑音生成用↓↓↓↓
#define sigmaN 0.5 				//ガウス雑音の分散σ^2

//確率密度関数用↓↓↓↓
#define range 1000				//試行範囲
#define w 0.05					//範囲の幅


//BER特性用↓↓↓↓
#define c 2 					//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1   				

////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void baseband(const double time[] , double complex base[]);
void boxmuller(double complex noise[],double mean,double variance);


int main(void)
{	
	int i,j; 

///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 
	
//////////////////////// 時間波形，スペクトル表示用の時間と周波数の設定　　/////////////////////
	double t[N]; 						//時間t
	double f[N];						//周波数f
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}

	
////////////////////　　ガウス雑音の生成　　////////////////////////
	double complex *awgn;				//ガウス雑音
	double complex *fft_awgn;			//ガウス雑音のスペクトル
	
	awgn=(double complex *)malloc(sizeof(double complex)*N);		//マロック関数により領域拡張
	fft_awgn=(double complex *)malloc(sizeof(double complex)*N);	//マロック関数により領域拡張
		
	
	
	boxmuller(awgn,0.0,sigmaN);		//ガウス雑音の生成
	FFT(awgn,fft_awgn);				//ガウス雑音のスペクトル算出　
	
////////////////////　　フェージングの作成　　///////////////////
	double complex *Rfading;			//レイリーフェージングの時間波形
	double complex *fft_Rfading;		//レイリーフェージングのスペクトル
	double complex *F_tmp;				//(レイリーフェージングのスペクトル)×(ガウス雑音のスペクトル)
	double *Pf_Rfading;					//レイリーフェージングの電力スペクトル密度
	double sigmaR=1.0;					//σ_R^2　電力スペクトル密度の総電力

	////////////// マロック関数による領域拡張　//////////////////
	Rfading=(double complex *)malloc(sizeof(double complex)*N);
	fft_Rfading=(double complex *)malloc(sizeof(double complex)*N);
	F_tmp=(double complex *)malloc(sizeof(double complex)*N);
	Pf_Rfading=(double *)malloc(sizeof(double)*N);
	
	/////////////　電力スペクトル密度　///////////////
	
	double Pt_ray=0.0;		//時間領域の総電力
	double Pf_ray=0.0;		//周波数領域の総電力
	double Pf_awgn=0.0;		//ガウス雑音の周波数領域での総電力
	double PP_ray=0.0;		//電力スペクトル密度の総和
//f-fcの部分について→f[i]=(i-N/2)delta_f=(i-N/2)(1/Ts*K)=
	for(i=0;i<N;i++){
		if(fd>fabs(f[i]))	//sqrt()の中身が負にならないように
			Pf_Rfading[i]=sigmaR/M_PI/(sqrt(fd*fd-f[i]*f[i]));	
		else 
			Pf_Rfading[i]=0.0;
		
		PP_ray+=Pf_Rfading[i];		//総電力=σR^2
	}
	
	printf("Freq\tF_Amp\n");
	for(i=0;i<N;i++){
		Pf_Rfading[i]=Pf_Rfading[i]/(PP_ray/N);					//電力スペクトル密度の平均電力を1Wに正規化している
		fft_Rfading[i]=sqrt(Pf_Rfading[i]);						//平均電力が1Wになった電力スペクトル密度のルートをとる
		F_tmp[i]=fft_awgn[i]*fft_Rfading[i];					//レイリーフェージングのスペクトル×AWGNのスペクトル
		Pf_ray+=cabs(fft_Rfading[i])*cabs(fft_Rfading[i]);		//レイリーフェージングの振幅スペクトルの平均電力 1Wになっている．112行目で確かめ
		Pf_awgn+=cabs(fft_awgn[i])*cabs(fft_awgn[i]);			//ガウス雑音の周波数領域の電力
		// printf("%f\t%e\n",f[i],cabs(fft_Rfading[i]));
	}

	///////　(レイリーフェージングのスペクトル)×(ガウス雑音のスペクトル)を時間領域に戻す  ///////
	IFFT(F_tmp,Rfading);
	
	for(i=0;i<N;i++)
		Pt_ray+=cabs(Rfading[i])*cabs(Rfading[i]);
	
	
	///////////////////　　電力の確認しよう！！　　//////////////////////////////////
	// printf("PP_ray\tPf_awgn/N\tPf_ray/N\tPt_ray/N\t\n");
	// printf("%f\t%f\t%f\t%f\n",PP_ray,Pf_awgn/N,Pf_ray/N,Pt_ray/N);		//平均電力の確かめ
	
	// putchar('\n');
	/* printf(" 時間　　　　　振幅波形　　　　位相波形    周波数　　　　　振幅スペクトル\n");
	for(i=0;i<N;i++){
		printf("%f\t%e\t%e\t%f\t%e\n",t[i],cabs(Rfading[i]),atan2(creal(Rfading[i]),cimag(Rfading[i]))/M_PI,f[i],cabs(fft_Rfading[i]));
	} */
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



//　　　　Box-Muller法によるガウス雑音生成　　　　　　//
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