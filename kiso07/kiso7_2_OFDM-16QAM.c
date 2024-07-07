////// OFDM-16QAMのAWGN通信路におけるBER特性 ///////////////

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1						//波数
#define s 1						//オーバーサンプル比
#define K 256					//
#define N (s*K)					//全サンプル数
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1シンボル時間
#define fc a/Ts 				//搬送波周波数fc
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長
//BER特性用↓↓↓↓
#define bit 4.0 				//1シンボル当たりのビット数 QPSKの場合は2
#define BnTs 1   				
#define loop pow(10,4)			//試行回数

////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void QAM16(double complex QAM[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void SP(double complex Serial[],double complex Parallel[]);
void PS(double complex Parallel[],double complex Serial[]);

int main(void)
{	
	int i,j,m; 
	
///////////////////////      乱数の種を設定     ///////////////
	// srand((unsigned int)time(NULL)); 

	
////////////////////// 時間と周波数の設定　　//////////////////////
	double t[N]; 						//時間t    0～K*Tsまでの時間
	double f[N];						//周波数f  -Ts/s～Ts/sまでの周波数
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}
	
	double complex sF1_16QAM[N];	//16QAM信号
	double complex sF2_16QAM[N];	//先頭サンプルのみ16QAM信号
	double complex sF[N];			//S/P変換後の信号
	double complex st[N];			//IFFT後の信号	
	double complex awgn[N];			//ガウス雑音
	double complex rt[N];			//
	double complex rF[N];			//
	double complex ro[N];			//
	double amp;						//16QAM信号の振幅
	
	amp=1/sqrt(10.0*s);
	
	QAM16(sF1_16QAM);				//16QAM信号生成
	inter(sF1_16QAM,sF2_16QAM);		//K本のサブキャリア=先頭サンプルのみ抜き出す
	SP(sF2_16QAM,sF);
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
	for(EbN0=0; EbN0<=15; EbN0++){
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
				printf("%f\t%f\t%f\n",t[i],creal(st[i]),creal(rt[i]));
			
			putchar('\n');
			printf("Freq\t|s(f)|\t|r(f)|\n");
			for(i=0;i<N;i++)
				printf */("%f\t%e\t%e\n",f[i],cabs(sF[i]),cabs(rF[i]));
			
			// putchar('\n');
			// for(i=0;i<N;i++)
				// printf("%f\t%f\t%f\t%e\n",t[i],creal(ro[i]),f[i],cabs(ro[i]));
			
			
			
			//// 誤り判定(先頭点のみ)  ////////////////////////////////////////////////////////////////////////
			for(j=0;j<K;j++){
			///////////////////////////  I軸判定　　//////////////////////////////////////////
				////////////    (00)伝送時の判定   /////////////
				if(creal(sF2_16QAM[j*s])>=2*amp){
					if(0<=creal(ro[j*s]) && creal(ro[j*s])<2*amp){	
						ci++;	//(00)伝送時(01)と誤る
					}else if(-2*amp<=creal(ro[j*s]) && creal(ro[j*s])<0){
						ci+=2;	//(00)伝送時(11)と誤る
					}else if(creal(ro[j*s])<=-2*amp){
						ci++;	//(00)伝送時(10)と誤る
					}
				}
					
				////////////    (01)伝送時の判定   /////////////
				if(0<=creal(sF2_16QAM[j*s]) && creal(sF2_16QAM[j*s])<2*amp){
					if(2*amp<=creal(ro[j*s])){	
						ci++;	//(01)伝送時(00)と誤る
					}else if(-2*amp<=creal(ro[j*s]) && creal(ro[j*s])<0){
						ci++;	//(01)伝送時(11)と誤る
					}else if(creal(ro[j*s])<=-2*amp){
						ci+=2;	//(01)伝送時(10)と誤る
					}	
				}
					
				////////////    (11)伝送時の判定   /////////////
				if(-2*amp<=creal(sF2_16QAM[j*s]) && creal(sF2_16QAM[j*s])<0){
					if(2*amp<=creal(ro[j*s])){	
						ci+=2;	//(11)伝送時(00)と誤る
					}else if(0<=creal(ro[j*s]) && creal(ro[j*s])<2*amp){
						ci++;	//(11)伝送時(01)と誤る
					}else if(creal(ro[j*s])<=-2*amp){
						ci++;	//(11)伝送時(10)と誤る
					}	
				}
					
				////////////    (10)伝送時の判定   /////////////
				if(creal(sF2_16QAM[j*s])<=-2*amp){
					if(2*amp<=creal(ro[j*s])){	
						ci++;	//(10)伝送時(00)と誤る
					}else if(0<=creal(ro[j*s]) && creal(ro[j*s])<2*amp){
						ci+=2;	//(10)伝送時(01)と誤る
					}else if(-2*amp<=creal(ro[j*s]) && creal(ro[j*s])<0){
						ci++;	//(10)伝送時(11)と誤る
					}	
				}
					
			///////////////////////////  Q軸判定　　//////////////////////////////////////////
				////////////    (00)伝送時の判定   /////////////
				if(cimag(sF2_16QAM[j*s])>=2*amp){
					if(0<=cimag(ro[j*s]) && cimag(ro[j*s])<2*amp){	
						cq++;	//(00)伝送時(01)と誤る
					}else if(-2*amp<=cimag(ro[j*s]) && cimag(ro[j*s])<0){
						cq+=2;	//(00)伝送時(11)と誤る
					}else if(cimag(ro[j*s])<=-2*amp){
						cq++;	//(00)伝送時(10)と誤る
					}
				}
					
				////////////    (01)伝送時の判定   /////////////
				if(0<=cimag(sF2_16QAM[j*s]) && cimag(sF2_16QAM[j*s])<2*amp){
					if(2*amp<=cimag(ro[j*s])){	
						cq++;	//(01)伝送時(00)と誤る
					}else if(-2*amp<=cimag(ro[j*s]) && cimag(ro[j*s])<0){
						cq++;	//(01)伝送時(11)と誤る
					}else if(cimag(ro[j*s])<=-2*amp){
						cq+=2;	//(01)伝送時(10)と誤る
					}	
				}
					
				////////////    (11)伝送時の判定   /////////////
				if(-2*amp<=cimag(sF2_16QAM[j*s]) && cimag(sF2_16QAM[j*s])<0){
					if(2*amp<=cimag(ro[j*s])){	
						cq+=2;	//(11)伝送時(00)と誤る
					}else if(0<=cimag(ro[j*s]) && cimag(ro[j*s])<2*amp){
						cq++;	//(11)伝送時(01)と誤る
					}else if(cimag(ro[j*s])<=-2*amp){
						cq++;	//(11)伝送時(10)と誤る
					}	
				}
					
				////////////    (10)伝送時の判定   /////////////
				if(cimag(sF2_16QAM[j*s])<=-2*amp){
					if(2*amp<=cimag(ro[j*s])){	
						cq++;	//(10)伝送時(00)と誤る
					}else if(0<=cimag(ro[j*s]) && cimag(ro[j*s])<2*amp){
						cq+=2;	//(10)伝送時(01)と誤る
					}else if(-2*amp<=cimag(ro[j*s]) && cimag(ro[j*s])<0){
						cq++;	//(10)伝送時(11)と誤る
					}	
				}
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

/////////////////////////////  16QAM信号の生成　　///////////////////////////
void QAM16(double complex QAM[]){
	int i,j;
	double b1,b2,b3,b4;		//乱数用変数
	double amp;				//振幅A
	double s1I,s1Q;			//I軸，Q軸の振幅
	
	////////////// 16QAM信号生成   ////////////////////
	amp=1/sqrt(10.0*s);		//振幅Aを算出
	
	for(i=0;i<K;i++){
		b1=rand()%2;	b2=rand()%2;		//I軸の乱数
		b3=rand()%2;	b4=rand()%2;		//Q軸の乱数
		
		////////　　同相成分決定   ////////////
		if(b1==0 && b2==0) s1I=3*amp;		//(00)なら同相成分は3A
		if(b1==0 && b2==1) s1I=amp;			//(01)なら同相成分はA
		if(b1==1 && b2==1) s1I=-amp;		//(11)なら同相成分は-A
		if(b1==1 && b2==0) s1I=-3*amp;		//(10)なら同相成分は-3A
		
		////////　　直交成分決定   ////////////
		if(b3==0 && b4==0) s1Q=3*amp;		//(00)なら直交成分は3A
		if(b3==0 && b4==1) s1Q=amp;			//(01)なら直交成分はA
		if(b3==1 && b4==1) s1Q=-amp;		//(11)なら直交成分は-A
		if(b3==1 && b4==0) s1Q=-3*amp;		//(10)なら直交成分は-3A
		
		for(j=0;j<s;j++)
			QAM[i*s+j]=s1I+I*s1Q;			//1シンボルごとにOS比分コピー
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












