#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//乱数の時間変化のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム

#define Ts 0.001				//1シンボル時間
#define a 3						//波数
#define s 2						//オーバーサンプル比
#define symble 256				//シンボル数
#define f0 (a/Ts)				//変調周波数
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define N (s*symble)			//全サンプル数
#define delta_f (1/(Ts*symble))	//Δf=1/T Δf:スペクトルの幅  T(=Ts*symble):時間長
	
//離散フーリエ変換(Discrete Fourier Transform)を行う関数//
void dft(double complex signal[] , double complex spectrum[]){
	for(int k=0; k<=N-1; k++){
		spectrum[k]=0;
	}
	
	for(int k=0; k<=N-1; k++){
		for(int n=0; n<=N-1; n++){
			spectrum[k]+=(1/sqrt((double)N))*signal[n]*cexp(-I*2*M_PI*(k-N/2)*(n-N/2)/N);
		//　　本当はkやnに-N/2～を入れたいが，C言語の配列は0からなのでk-N/2やn-N/2にする　　//
		}
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

int main(void)
{
	int i,j;  
	double I_t[N];						//複素ベースバンド信号の同相成分
	double Q_t[N];						//複素ベースバンド信号の直交成分
	double complex u_t[N];				//複素ベースバンド信号
	double complex s_t[N];				//搬送帯域信号
	double t=0.0; 						//時間t
	
	//      ランダム関数初期化      //
	srand((unsigned int)time(NULL)); 
	
	//　　　　　　複素ベースバンド信号と搬送帯域信号の生成　　　　　//
	double tmp_I;	//同相成分Iを一時的に格納
	double tmp_Q;	//直交成分Qを一時的に格納
	
	for(i=0; i<symble; i++){
		tmp_I=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
		tmp_Q=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
		
		for(j=0; j<s; j++){
			I_t[s*i+j]=tmp_I;				//同相成分I(t)にs個(サンプル数)だけ格納
			Q_t[s*i+j]=tmp_Q;				//直交成分Q(t)にs個(サンプル数)だけ格納
			u_t[s*i+j]=tmp_I+tmp_Q*I;		//複素ベースバンド信号u(t)=I(t)+jQ(t)
			
			s_t[s*i+j]=creal(u_t[s*i+j])*cos(2*M_PI*f0*t)-cimag(u_t[s*i+j])*sin(2*M_PI*f0*t);	//搬送帯域信号s(t)
			
			/*
			printf("%f  ", t);
			printf("%f  ", creal(u_t[s*i+j]));
			printf("%f  ", cimag(u_t[s*i+j]));
			printf("%f  \n", s_t[s*i+j]);
			*/
			
			t=t+T_sample;
		}
	}
	
	
	//		離散フーリエ変換		//
	double complex spectrum_base[N];	//ベースバンド信号のスペクトル
	double complex spectrum_carrier[N];	//搬送帯域信号のスペクトル
	
	dft(u_t,spectrum_base);
	dft(s_t,spectrum_carrier);
	
	printf("\n\n");
	
	
	// for(int i=0; i<N; i++){
	// printf("%f   " , creal(spectrum_base[i]));
	// printf("%f   " , cimag(spectrum_base[i]));
	// printf("%f\n" , cabs(spectrum_base[i]));
	// }
	
	// printf("\n");
	
	// for(int i=0; i<N; i++){
	// printf("%f   " , creal(spectrum_carrier[i]));
	// printf("%f   " , cimag(spectrum_carrier[i]));
	// printf("%f\n" , cabs(spectrum_carrier[i]));
	// }	
	
	//		周波数軸に変換 (spectrum[0]を周波数-N/2*Δfに割り当て)		//
	double f;			//周波数
	
	printf("\n\n");
	printf("離散フーリエ変換によるスペクトルの表示\n");
	printf("周波数　   　　ベース振幅　       　　搬送波スペクトル \n");
	for(int i=0; i<N; i++){
		f=(i-N/2)*delta_f;
		//  周波数   //
		printf("%.1f      " , f);
		
		//　　ベースバンド信号のスペクトル　　//
		// printf("%.15f   " , creal(spectrum_base[i]));
		// printf("%.15f   " , cimag(spectrum_base[i]));
		printf("%.15f       " , cabs(spectrum_base[i]));
		
		//　　　搬送帯域信号のスペクトル　　　//
		// printf("%.15f   " , creal(spectrum_carrier[i]));
		// printf("%.15f   " , cimag(spectrum_carrier[i]));
		printf("%.15f\n" , cabs(spectrum_carrier[i]));

	}
	
	//　               　　電力の計算（時間領域，周波数領域） 　                 　//
	double Pt_base=0.0;					//時間領域のベースバンド信号の電力
	double Pt_carrier=0.0;				//時間領域の搬送帯域信号の電力
	double Pf_base=0.0;					//周波数領域のベースバンド信号の電力
	double Pf_carrier=0.0;				//周波数領域の搬送帯域信号の電力
	
	printf("\n\n");
	printf("時間領域と周波数領域での電力の表示 \n");
	for(int i=0; i<N; i++){
			Pt_base+=cabs(u_t[i])*cabs(u_t[i]);									//ベースバンド信号の電力（時間領域）
			Pt_carrier+=cabs(s_t[i])*cabs(s_t[i]);								//搬送帯域信号の電力（時間領域）
			Pf_base+=cabs(spectrum_base[i])*cabs(spectrum_base[i]);				//ベースバンド信号の電力（周波数領域）
			Pf_carrier+=cabs(spectrum_carrier[i])*cabs(spectrum_carrier[i]);	//搬送帯域信号の電力（周波数領域）
	}
	
	printf("Pt_base=%.2f [W]  Pt_carrier=%.2f [W] \n",Pt_base,Pt_carrier);
	printf("Pf_base=%.2f [W]  Pf_carrier=%.2f [W] \n",Pf_base,Pf_carrier);
	
	
	
	
	
	//	　離散フーリエ逆変換　　//
	
	double complex tmp_base[N];
	double complex tmp_carrier[N];
	
	idft(spectrum_base,tmp_base);
	idft(spectrum_carrier,tmp_carrier);
	
	t=0.0;
	//            本当に信号が戻っているかの確認　　               //
	printf("\n\n");
	printf("ベースバンド信号が離散フーリエ逆変換で元に戻っているかの確認\n");
	printf("  時間t    u(t)のI軸   tmpのI軸   u(t)のQ軸     tmpのQ軸\n");
	for(int i=0; i<N; i++){
		printf("%f   " , t);
		printf("%f   " , creal(u_t[i]));
		printf("%f   " , creal(tmp_base[i]));
		printf("%f   " , cimag(u_t[i]));
		printf("%f\n" , cimag(tmp_base[i]));
		t+=T_sample;
	}
	t=0.0;
	printf("\n\n");
	printf("搬送帯域信号が離散フーリエ逆変換で元に戻っているかの確認\n");
	printf("  時間t     s(t)のI軸   tmpのI軸   s(t)のQ軸     tmpのQ軸\n");
	for(int i=0; i<N; i++){
		printf("%f   " , t);
		printf("%f   " , creal(s_t[i]));
		printf("%f   " , creal(tmp_carrier[i]));
		printf("%f   " , cimag(s_t[i]));
		printf("%f\n" , cimag(tmp_carrier[i]));
		t+=T_sample;
	}
	
	
	return 0;	
}