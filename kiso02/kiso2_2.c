#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//乱数の時間変化のプログラム,clock関数使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define Ts 0.001				//1シンボル時間
#define a 1						//波数
#define s 1						//オーバーサンプル比
#define symble 256				//シンボル数
#define f0 (a/Ts)				//変調周波数
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define N (s*symble)			//全サンプル数
//sやsymbleの値を変えたら，bit_lengthの値も変える  //
#define bit_length 12			//ビット長 N=2^bit_length
#define delta_f (1/(Ts*symble))	//Δf=1/T Δf:スペクトルの幅  T(=Ts*symble):時間長


//離散フーリエ変換(Discrete Fourier Transform)を行う関数//
void dft(double complex signal[] , double complex spectrum[]){
	int k, n;
	for(k=0; k<=N-1; k++){
		spectrum[k]=0;
	}
	
	for(k=0; k<=N-1; k++){
		for(n=0; n<=N-1; n++){
			spectrum[k]+=(1/sqrt((double)N))*signal[n]*cexp(-I*2*M_PI*(k-N/2)*(n-N/2)/N);
		//　　本当はkやnに-N/2～を入れたいが，C言語の配列は0からなのでk-N/2やn-N/2にする　　//
		}
	}
}

// //　　　　入力を並び替えコピーする関数　　　//
// void bit_reverse(int bit[]){
	
// 	bit[0] = 0;
// 	int i,j;
// 	int k = 1;
// 	for(i = 0; i < bit_length; i++){
// 		for(j = 0; j < k; j++){
// 			bit[j] *= 2;
// 		}
// 		for(j = 0; j < k; j++){
// 			bit[j+k] = bit[j] + 1;
// 		}
// 		k *= 2;
// 	}
// }

// //　　　　数値を入れ替える関数　　　　//
// 				// sample                tmp
// void swap(double complex x[] , double complex y[]){
// 	int i;
// 	int input[N];				//入力並び替えのための数字をもらう

// 	bit_reverse(input);			//ビット反転した数値を入れる
	
// 	for(i=0; i<N; i++)
// 		y[i]=x[input[i]];	
// }

// //　高速フーリエ変換(Fast Fourier Transform) //
// void fft(double complex sample[] , double complex spectrum[]){
// 	int k,n;
	
// 	double complex tmp[N];
// 	double complex X[N];
// 	int count=0;		//配列の要素用カウンタ
// 	int m;				//角度調整用カウンタ
	
// 	// bit_reverse関数とswap関数に入れる  tmp[N]に並び替えられた入力が入る
// 	swap(sample,tmp);
	
// 	for(k=2; k<=N; k*=2){	   //N=2^nとしたときn回だけ繰り返す．(バタフライがn段ある)				
// 		for(n=0; n<N; n+=k){	
// 			m=0;

// 			//　　　バタフライ演算　　　//
// 			X[n]=tmp[n] + cexp(-I*2*M_PI*(m-N/2)/k) * tmp[n+(1<<count)];							//かっこいる？
// 			X[n+(1<<count)]=tmp[n] - cexp(-I*2*M_PI*(m-N/2)/k) * tmp[n+(1<<count)];	
// 			// X[n]=tmp[n]+cexp(-I*2*M_PI*m/k)*tmp[n+(1<<count)];							//かっこいる？
// 			// X[n+(1<<count)]=tmp[n]-cexp(-I*2*M_PI*m/k)*tmp[n+(1<<count)];	
			
// 			// 　コピー //
// 			tmp[n]=X[n];
// 			tmp[n+(1<<count)]=X[n+(1<<count)];
			
// 			//　　　初段以外はここのループに入る    //
// 			for(m=1; m<=(k/2-1); m++){
// 			X[n+m]=tmp[n+m]+cexp(-I*2*M_PI*(m-N/2)/k)*tmp[n+m+(1<<count)];							//かっこいる？
// 			X[n+m+(1<<count)]=tmp[n+m]-cexp(-I*2*M_PI*(m-N/2)/k)*tmp[n+m+(1<<count)];
// 			// X[n+m]=tmp[n+m]+cexp(-I*2*M_PI*m/k)*tmp[n+m+(1<<count)];							//かっこいる？
// 			// X[n+m+(1<<count)]=tmp[n+m]-cexp(-I*2*M_PI*m/k)*tmp[n+m+(1<<count)];
			
// 			tmp[n+m]=X[n+m];
// 			tmp[n+m+(1<<count)]=X[n+m+(1<<count)];
// 				}
// 		}	
// 		count++;
// 	}
	
// 	for(k=0; k<=N-1; k++){
// 		spectrum[k]=sqrt(1/(double)N)*X[k];
// 	}
// }

//　　　　先輩のFFT　　　　　　//
void FFT(double complex Time[] , double complex Freq[]){

	int i,k,m,j;
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
		// 値の写し
		for(i=0; i<N; i++){
			X[i]=Y[i];
		}
		// Yのリセット
		for(i=0; i<N; i++)
			Y[i]=0.0+0.0*I;
	}
	
	for(i=0;i<N;i++)
		Freq[i]=X[i]/sqrt((double)N);
}

int main(void)
{
	int i,j;  
	double I_t[N];						//複素ベースバンド信号の同相成分
	double Q_t[N];						//複素ベースバンド信号の直交成分
	double complex u_t[N];				//複素ベースバンド信号
	complex s_t[N];				//搬送帯域信号
	double t=0.0; 						//時間t
	double f;							//周波数f
	long cpu_time;						//ストップウォッチ(処理時間を計測)
	
	//      ランダム関数初期化      //
	srand((unsigned int)time(NULL)); 
	
	//////////////　　　　　　複素ベースバンド信号と搬送帯域信号の生成  　　　////////////////
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
			
			
			printf("%f  ", t);
			printf("%f  ", creal(u_t[s*i+j]));
			printf("%f  ", cimag(u_t[s*i+j]));
			printf("%f  \n", s_t[s*i+j]);
			
			t=t+T_sample;
		}
	}
	
	////////////////////      離散フーリエ変換        //////////////////////
	double complex spectrum_base[N];	//ベースバンド信号のスペクトル
	double complex spectrum_carrier[N];	//搬送帯域信号のスペクトル
	
	dft(u_t,spectrum_base);
	dft(s_t,spectrum_carrier);
	
	printf("\n\n");
	printf("      離散フーリエ変換によるスペクトルの表示\n");
	printf("周波数　   　　ベース振幅　       　　搬送波スペクトル \n");
	
	for(i=0; i<N; i++){
		f=(i-N/2)*delta_f;
		//  周波数   //
		printf("%.1f      " , f);
		
		//　　ベースバンド信号のスペクトル　　//
		// printf("%.5f   " , creal(spectrum_base[i]));
		// printf("%.5f   " , cimag(spectrum_base[i]));
		printf("%.15f       " , cabs(spectrum_base[i]));
		
		//　　　搬送帯域信号のスペクトル　　　//
		// printf("%.5f   " , creal(spectrum_carrier[i]));
		// printf("%.5f   " , cimag(spectrum_carrier[i]));
		printf("%.15f\n" , cabs(spectrum_carrier[i]));
	}
	

	/////////////////////　     　　高速フーリエ変換　　　   　    //////////////////
	double complex fft_base[N];
	double complex fft_carrier[N];
	
	FFT(u_t,fft_base);
	FFT(s_t,fft_carrier);
	
	// printf("\n\n");
	printf("\n\n");
	printf("     高速フーリエ変換によるスペクトルの表示\n");
	printf("周波数　   　　ベース振幅　       　　搬送波スペクトル \n");
	for(i=0; i<N; i++){
		f=(i-N/2)*delta_f;
		//  周波数   //
		printf("%.1f      " , f);
		
		//　　ベースバンド信号のスペクトル　　//
		// printf("%.5f   " , creal(fft_base[i]));
		// printf("%.5f   " , cimag(fft_base[i]));
		printf("%.15f       " , cabs(fft_base[i]));
		
		//　　　搬送帯域信号のスペクトル　　　//
		// printf("%.5f   " , creal(fft_carrier[i]));
		// printf("%.5f   " , cimag(fft_carrier[i]));
		printf("%.15f\n" , cabs(fft_carrier[i]));
	}
	
	/////////////     処理にかかったCPU時間　　　　　　//////////////
	cpu_time=clock();
	
	putchar('\n');
	printf("CPU時間 %f秒 \n" , (double)cpu_time/CLOCKS_PER_SEC);
	
	return 0;
}