///////  マルチパスフェージングの作成（OFDM伝送はkiso8_2.cで作成）  //////////////

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1						//波数
#define s 4						//オーバーサンプル比
#define K 64					//サブキャリア数
#define N (s*K)					//全サンプル数
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1シンボル時間
#define fc a/Ts 				//搬送波周波数fc
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

//OFDM特有のパラメータ↓↓↓↓
#define L 4						//L波マルチパス
#define dB 1					//何dB減衰
#define GI 8					//ガードインターバルのサンプル数


//BER特性用↓↓↓↓
#define bit 2.0 				//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1   				
#define loop 1			//試行回数

////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void QPSK(double complex base[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void SP(double complex Serial[],double complex Parallel[]);
void PS(double complex Parallel[],double complex Serial[]);
void qs_Rayleigh(double complex qsRfading[]);

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
	
///////////////////// 準静的dB減衰L波マルチパスフェージングモデルの作成  ////////////////////////
	double complex ht_multifd[N];		//マルチパスフェージングのインパルス応答
	double complex Hf_multifd[N];		//マルチパスフェージングの伝達関数
	double P_multifd[L];				//各パスの電力
	double complex awgn1[1];			
	double Psum1=0.0;						//全パスの電力の総和
	double Psum2=0.0;
	double Psum3=0.0;			//正規化H(f)の総電力
	//マルチパスフェージング作成
	
	P_multifd[0]=(1-pow(10,-(dB/10.0)))/(1-pow(10,-(dB*L/10.0)));		//各パスのインパルス応答の電力の総和が1Wとなるためのh0の値
	/* for(i=0;i<K;i++){		
		if(i<L){		//L波分作成
			P_multifd[i]=P_multifd[0]*pow(10,-(dB*i/10.0));		//電力を減衰させてゆく P_multifd[0]=P_multifd[0]*1となる．
			boxmuller(1,awgn1,0.0,P_multifd[i]/2.0);					//2σ^2=P_multifd[i]となるようにガウス雑音を作成	
			for(j=0;j<s;j++)	//s個同じ準静的フェージング作成
				ht_multifd[j+s*i]=awgn1[0];
		}else{
			for(j=0;j<s;j++)	
				ht_multifd[j+s*i]=0+0*I;
		}
	} */   //準静的レイリーフェージングをs個同じawgnを入れたが，Ts毎のみでよかった
	
	
	for(i=0;i<K;i++){
		for(j=0;j<s;j++){
		if(i<L && j==0){		//L波分作成
			P_multifd[i]=P_multifd[0]*pow(10,-(dB*i/10.0));		//電力をYdBずつ減衰させてゆく P_multifd[0]=P_multifd[0]*1となる．
			boxmuller(1,awgn1,0.0,P_multifd[i]/2.0);			//2σ^2=P_multifd[i]となるようにガウス雑音を作成	
			ht_multifd[j+s*i]=awgn1[0];
		}else
			ht_multifd[j+s*i]=0+0*I;
		}
	}
	
	//確認用
	for(i=0;i<N;i++)
		printf("%f\t%f\n",t[i],cabs(ht_multifd[i]));

	FFT(ht_multifd,Hf_multifd);		//マルチパスフェージングの伝達関数H(f)
	
	//正規化前のH(f)　　最初失敗したもの
	// for(i=0;i<N;i++)
		// printf("%f\t%f\n",f[i],cabs(Hf_multifd[i]));
	
	//P_multifdの総電力計算→1Wになってほしい
	for(i=0;i<L;i++)
		Psum2+=P_multifd[i];
	
	//H(f)の総電力計算→正規化の準備
	for(i=0;i<N;i++)
		Psum1+=cabs(Hf_multifd[i])*cabs(Hf_multifd[i]);
	
	//平均電力で正規化    　※総電力で正規化したら間違った振幅スペクトルが得られてしまった
	for(i=0;i<N;i++)
		Hf_multifd[i]=Hf_multifd[i]/sqrt(1.0/N);		//H(f)は電圧なので，平均電力のルートで正規化
	
	
	//H(f)の総電力計算
	for(i=0;i<N;i++)
		Psum3+=cabs(Hf_multifd[i])*cabs(Hf_multifd[i]);
	
	
	
	putchar('\n');
	printf("%dパス準静的レイリーフェージング，%ddB減衰モデル　s=%d\n",L,dB,s);
	for(i=0;i<N;i++)
		printf("%f\t%f\n",f[i],cabs(Hf_multifd[i]));
	
	printf("P_multifdの総和=%f W\t正規化後H(f)の平均電力=%f W\n",Psum2,Psum3/N);
	
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
/////////////////////////////// 準静的レイリーフェージングの生成　//////////////////////////////////////
void qs_Rayleigh(double complex qsRfading[]){
	int i;
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
