////////////////////準静的レイリーフェージング作成 ////////////////////////////////

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1						//波数
#define s 8						//オーバーサンプル比
#define K 16					//1フレーム長
#define N s*K					//全サンプル数
#define frame 100				//フレーム数

#define Ts 0.001				//1シンボル時間
#define fc a/Ts 				//搬送波周波数fc
#define fdTs 0.001				//fdTs
#define fd (fdTs/Ts)			//ドップラーシフト
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

//ガウス雑音生成用↓↓↓↓
#define sigmaN 0.5 				//ガウス雑音の分散σ^2

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
	
	//ファイルへの書き込み
	FILE *fp;
	char *fname="result_kiso5_3.csv";
	fp=fopen(fname, "w");
	if(fp==NULL){
	printf("%sファイルが開けません\n", fname);
		return (-1);
	}
		else{

		double complex *awgn;				//ガウス雑音
		double complex *Rfading;			//準静的レイリーフェージングの時間波形
		double *t2;							//準静的レイリーフェージング用の時間
		////////////// マロック関数による領域拡張　//////////////////
		awgn=(double complex *)malloc(sizeof(double complex)*N);		//配列の要素数をframeにしてもできる　m=0～frame-1までの計frame個だから			
		Rfading=(double complex *)malloc(sizeof(double complex)*(K*frame));	//しかしboxmuller関数がN個からなる配列用の関数なのでawgn[N]のままでやるN>frameが条件
		t2=(double *)malloc(sizeof(double)*(K*frame));
		
		////////////////////　　ガウス雑音の生成　　////////////////////////
		boxmuller(awgn,0.0,sigmaN);		//ガウス雑音の生成
			
		putchar('\n');
		printf(" 時間　　　　　振幅波形　　　　位相波形  \n");
		for(int m=0;m<frame;m++){		//100フレーム数繰り返す	
			////////////　　準静的レイリーフェージングの作成＝ガウス雑音を1フレームずつコピーする  ////////////////
			for(i=0;i<K;i++){	//1フレーム単位でコピー
				Rfading[i+m*K]=awgn[m];		//1フレーム，つまり16点分同じガウス雑音をコピーする，awgnはframe分配列があればOK
				t2[i+m*K]=Ts*(i+m*K);		//準静的フェージング1フレーム(16点)でK*Ts秒，つまり1点でTs秒割り振る
				fprintf(fp,"%f\t%e\t%e\n",t2[i+m*K],cabs(Rfading[i+m*K]),carg(Rfading[i+m*K])/M_PI);
			}
			///////////////  結果表示　　////////////////
			// for(i=0;i<K;i++)
				// printf("%f\t%e\t%e\n",t2[i+m*K],cabs(Rfading[i+m*K]),carg(Rfading[i+m*K])/M_PI);
		}
	////////　メモリを解放　/////////
	free(awgn);
	free(Rfading);
	free(t2);	
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