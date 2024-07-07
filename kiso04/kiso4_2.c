#include<stdio.h>
#include<stdlib.h> /*乱数*/
#include<time.h> /*乱数の時間変化*/
#include<math.h> /*三角関数の計算*/
#include<complex.h> /*複素数計算*/

#define K 256 /* シンボル数 */
#define s 8 /* OS比 */
#define Ts 0.001 /* 1シンボル時間 */
#define N (K*s) /* サンプル数 */
#define Tsample (Ts/s) /*1サンプル当たりの時間*/
#define delta_f (1/(Ts*K)) /*Δf=1/T Δf:スペクトルの幅 T(=Ts*SYMBOL):時間長*/
#define f0 1000 /* 帯域幅 */
#define alpha 0.5 /* ロールオフ率 */
#define EbN0 0 
#define bit 2 

/* 関数プロトタイプ宣言（最後に記述） */
void DFT(double complex signal[],double complex spectrum[]);
void IDFT(double complex spectrum[],double complex signal[]);
void QPSK(double complex BB[]);
void sample1(double complex BB[],double complex BB1[]);
void imp(double complex h[]);
void AWGN(double complex GZ[],double mean,double variance);
void tatami(double complex FBB[],double complex H[],double complex LPFBB[]);


int main(void){
	
	srand((unsigned)time(NULL)); /* 乱数初期化 */
	
	double complex h[N];
	double t[N];
	int n;
	double complex H[N]; /* 周波数応答 */
	double f[N]; /* 周波数軸 */
	double complex BB[N]; /* 複素ベースバンド信号 */
	double complex BB1[N]; /* 一点目だけ取り出した複素ベースバンド信号 */
	double complex FBB1[N];
	double complex LPFBB[N]; /* 複素ベースバンド信号と周波数応答の積 */
	double complex lpfBB[N]; /* LPF後の信号 */
	double complex GZ[N];
	double sigma2;
	double complex r0[N]; /* 受信信号（送信信号＋ガウス雑音） */
	double complex Fr01[N];
	double complex LPFr0[N];
	double complex r[N];
	double complex FBB[N];
	double complex Fr0[N];	
	
	imp(h);

	for(n=0;n<N;n++){
		
		t[n]=Tsample*(n-N/2); /* 時間軸作る */

	// printf("%f\t%f\n",t[n],h[n]); /* インパルス応答出力 */
	}
	
	DFT(h,H);
	
	for(n=0;n<N;n++){
		
		f[n]=delta_f*(n-N/2); /* 周波数軸作る */
		
		printf("%f\t%f\t%f\n",f[n],creal(H[n]),cimag(H[n])); /* 周波数応答出力 */
		
	}
	
	for(n=0;n<N;n++){
		
		H[n]=sqrt(cabs(H[n])); /* フィルタをルートにする */
		
		// printf("%f\t%f\n",f[n],H[n]);
	}
	
	IDFT(H,h);
	
	// for(n=0;n<N;n++){
	// printf("%f\t%f\n",t[n],h[n]); /* インパルス応答出力 */
	// }
	/* 複素ベースバンド信号をフィルタリング */
	
	
	QPSK(BB);

	// printf("複素ベースバンド信号\n");
	// for(n=0;n<N;n++){
		// printf("%f\t%f\t%f\n",t[n],creal(BB[n]),cimag(BB[n])); /* 複素ベースバンド信号の出力 */
	// }
	

	sample1(BB,BB1);
		
	// for(n=0;n<N;n++){
		// printf("%f\t%f\t%f\n",t[n],creal(BB1[n]),cimag(BB1[n]));
	// }
	

	DFT(BB1,FBB1);
	tatami(FBB1,H,LPFBB);
	IDFT(LPFBB,lpfBB);

	//雑音生成
	sigma2=1.0/(2*s*bit)*pow(10.0,-1.0*EbN0/10); /* 分散を計算 */
	AWGN(GZ,0,sigma2);
	
	//雑音を信号に足す
	for(n=0;n<N;n++){
		r0[n]=lpfBB[n]+GZ[n];
		
		// printf("%f\t%f\t%f\n",t[n],creal(r0[n]),cimag(r0[n]));
	}
	
	// double complex r01[N];
	
	// sample1(r0,r01);
	

	DFT(r0,Fr01);
	tatami(Fr01,H,LPFr0);
	IDFT(LPFr0,r);
	
	/* s0を出力 */
	// printf("s0の波形\n");
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\t%f\n",t[n],creal(BB[n]),cimag(BB[n]));
		
	// }
	
	// /* sを出力 */
	// printf("sの波形\n");
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\t%f\n",t[n],creal(lpfBB[n]),cimag(lpfBB[n]));
		
	// }
	
	// /* r0を出力 */
	// printf("r0の波形\n");
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\t%f\n",t[n],creal(r0[n]),cimag(r0[n]));
		
	// }
	
	// /* rを出力 */
	// printf("rの波形\n");
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\t%f\n",t[n],creal(r[n]),cimag(r[n]));
		
	// }
	
	
	/* スペクトル */
	
	
	DFT(BB,FBB);
	
	// /* s0を出力 */
	// printf("s0のスペクトル\n");
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\n",f[n],cabs(FBB[n]));
		
	// }
	
	// /* sを出力 */
	// printf("sのスペクトル\n");
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\n",f[n],cabs(LPFBB[n]));
		
	// }
	
	DFT(r0,Fr0);
	// /* r0を出力 */
	// printf("r0のスペクトル\n");
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\n",f[n],cabs(Fr0[n]));
		
	// }
	
	// /* rを出力 */
	// printf("rのスペクトル\n");
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\n",f[n],cabs(LPFr0[n]));
		
	// }
	
	/* I軸 */
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\t%f\t%f\t%f\n",t[n],creal(BB[n]),creal(lpfBB[n]),creal(r0[n]),creal(r[n]));
		
	// }
	
	// /* Q軸 */
	// for(n=0;n<N;n++){
		
		// printf("%f\t%f\t%f\t%f\t%f\n",t[n],cimag(BB[n]),cimag(lpfBB[n]),cimag(r0[n]),cimag(r[n]));
		
	// }
	
	// /* スペクトル */
	// for(n=0;n<N;n++){
		
		// printf("%f\t%.15f\t%.15f\t%.15f\t%.15f\n",f[n],cabs(FBB[n]),cabs(LPFBB[n]),cabs(Fr0[n]),cabs(LPFr0[n]));
		
	// }
	
	
	return 0;
}


// 離散フーリエ変換(DFT) 標本化された時間信号から標本化されたスペクトル信号に
void DFT(double complex signal[],double complex spectrum[]){
	
	int k,n; //kはスペクトル信号，nは時間信号のカウンタ
	
	// spectrum[]にN個の初期値を入れる
	for(k=0;k<=N-1;k++){
		spectrum[k]=0;
	}

	for(k=0;k<=N-1;k++){
		//スペクトル信号の離散値のそれぞれの値を求める(すべてのkを計算)
		for(n=0;n<=N-1;n++){
			//nを0~N-1まで変えた時の和
			spectrum[k]+=(1/sqrt((double) N))*signal[n]*cexp((-I*2*M_PI*(k-N/2)*(n-N/2))/N);
			//k-N/2,n-N/2にして-N/2~N/2にする
		}
	}
}

// 離散フーリエ逆変換(IDFT)
void IDFT(double complex spectrum[],double complex signal[]){
	
	int k,n; //kはスペクトル信号，nは時間信号のカウンタ
	
	// signal[]にN個の初期値を入れる
	for(n=0;n<=N-1;n++){
		signal[n]=0;
	}

	for(n=0;n<=N-1;n++){
		//時間信号の離散値のそれぞれの値を求める
		for(k=0;k<=N-1;k++){
			//nを0~N-1まで変えた時の和
			signal[n]+=(1/sqrt((double) N))*spectrum[k]*cexp((I*2*M_PI*(k-N/2)*(n-N/2))/N);
			//k-N/2,n-N/2にして-N/2~N/2にする
		}
	}
}

// ベースバンド信号を作る
void QPSK(double complex BB[]){
	double A;
	int QPSKbit;
	int k;
	double It[N];
	double Qt[N];
	
	A=1/(sqrt(2*s));//1シンボルの電力を1にする。1サンプルの電力は1/s
	
	// I軸送信信号を作る
	for(k=0;k<K;k++){
		QPSKbit=-2*(rand()%2)+1;
		It[k]=A*QPSKbit;
	}
	
	// Q軸送信信号を作る
	for(k=0;k<K;k++){
		QPSKbit=-2*(rand()%2)+1;
		Qt[k]=A*QPSKbit;
	}
	
	int i;
	
	for(k=0;k<K;k++){
		for(i=0;i<s;i++){
			BB[s*k+i]=It[k]+I*Qt[k];
		}
	}
}


/*複素ベースバンド信号のそれぞれのシンボルの1点目だけを取り出す*/
void sample1(double complex BB[],double complex BB1[]){
	
	int k; /* シンボル数える */
	int i; /* 点数える */
	
	for(k=0;k<K;k++){
		BB1[k*s]=BB[k*s];
		for(i=1;i<s;i++){
			BB1[k*s+i]=0;
		}
	}
}

void imp(double complex h[]){
	int n;
	double t[N]; /* 時間 */
	
	for(n=0;n<N;n++){
		
		t[n]=Tsample*(n-N/2); /* インパルス応答の時間軸 */
		
	
		if(t[n]==0.0)
			h[n]=sqrt(N)*1.0; /* ルートN掛ける理由は謎 */
		
		else if(t[n]==Ts/(2.0*alpha)||t[n]==-Ts/(2.0*alpha))
			h[n]=sqrt(N)*M_PI/4.0*sin(M_PI*t[n]/Ts)/(M_PI*t[n]/Ts);
		
		else
			h[n]=sqrt(N)*sin(M_PI*t[n]/Ts)/(M_PI*t[n]/Ts)*cos(M_PI*alpha*t[n]/Ts)/(1.0-pow(2.0*alpha*t[n]/Ts,2));
	}
}

//雑音生成
void AWGN(double complex GZ[],double mean,double variance){

	int n; /* サンプルカウント */
		double u1,u2; /* (0,1]の一様乱数を入れます */
		double XI[N],XQ[N]; /* Box-Muller法の標準正規分布 */
		double k=2.0; /* ビット数 */
		
		
		
		/* 複素ガウス雑音(I軸とQ軸の時間波形を求める，振幅スペクトルを求める） */
		for(n=0;n<N;n++){
			
			/* (0,1]の一様乱数を作る */
			u1=(rand()+1.0)/(RAND_MAX+1.0); /* あとでlogに入れるから分母分子にそれぞれ1を足して0にならないようにする */
			u2=(rand()+1.0)/(RAND_MAX+1.0); /*int÷int=intになってしまうので片方doubleにしておく*/
			
			XI[n]=sqrt(-2.0*variance*log(u1))*cos(2*M_PI*u2)+mean;
			XQ[n]=sqrt(-2.0*variance*log(u1))*sin(2*M_PI*u2)+mean; /* 分散がσ^2の複素ガウス雑音 */
			
			GZ[n]=XI[n]+I*XQ[n]; /* 複素ガウス雑音 */
			
		}
}

//周波数領域の畳み込み
void tatami(double complex FBB[],double complex H[],double complex LPFBB[]){
	
	int n;
	
	for(n=0;n<N;n++){
		LPFBB[n]=FBB[n]*H[n];
	}
}	
	
