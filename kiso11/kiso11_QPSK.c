#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>

#define s 1	//OS比	Bn*Tsの値としても用いる
#define Ts 0.001//1シンボル時間
#define symbol 3200	//シンボル数
#define N (s*symbol)		//合計サンプリング数
#define BT 1
#define BIT 2
#define EbN0 12
#define data 4 //情報ビット長
#define code 7	//符号語長


/*================================== 関数プロトタイプ宣言 ===========================================*/
void DFT(double complex Freq[], double complex Time[]);
void IDFT(double complex Time[],double complex Freq[]);
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void QPSKbaseband(double complex base[]);
void SIXTQAMbaseband(double complex base[]);
void AWGN(double sigma,double complex cgn[]);

/*==============================================================================================*/

int main(void){
	
	/*------------------- 開始時間の表示 ----------------------------------------*/
	time_t timer;
    struct tm *date;
    char str[256];

    timer = time(NULL);          /* 経過時間を取得 */
    date = localtime(&timer);    /* 経過時間を時間を表す構造体 date に変換 */
    printf("%s\n", asctime(date));                    /* 構造体 date を文字列に変換して表示 */
    strftime(str, 255, "%Y, %B, %d, %A %p%I:%M:%S", date);
    //printf("%s\n\n\n", str);
	/*-----------------------------------------------------------------------------*/
	
	
	
	
//////////////////////////////////////////////////////////////////////////////////////////////	
	int i,j,ENdB,l;
	double sec,Hz;
		
	//情報ビットの発生// //QSPKベースバンド信号の発生//
	int Iaxis[symbol], Qaxis[symbol];
	double complex base[symbol], basedash[N];
	
	//AWGNの発生//
	double sigma;
	double complex gauss[N],basegauss[N];
	
	//復号//
	int Iaxdecode[symbol],Qaxdecode[symbol],Idecode[symbol],Qdecode[symbol];
	int syndrome[3];
	
	//誤り符合//
	double BER;
	int icount,qcount;
	
//////////////////////////////////////////////////////////////////////////////////////////////	
	
	/*--------------------- ファイル名 --------------------------------------*/
	printf("ファイル名:kiso11_QPSK.c\n\n\n");
	
	
	/*--------------------- パラメータの表示 ---------------------------------*/
	int *p,loop;
	p=&loop;	*p=10000;
	printf("s=%d\t", s);	printf("symbol=%d\t", symbol);	printf("loop=%d\n", loop);
	printf("===================================================================================\n");
	
	/*--------------------------------------------------------------------*/
	
	/*---------------------- 実行結果の出力表示時のラベル --------------------*/
	// printf("試行回数=%d\t",loop);
	// printf("Eb/N0[dB]\t s=%d\n", s);
	
	/*---------------------------------------------------------------*/
	
	
	
	srand((unsigned)time(NULL));
	
	
/*===========================================================================================================================*/
	for(ENdB=0;ENdB<EbN0;ENdB++){
		icount=0;
		qcount=0;
		
		for(l=0;l<loop;l++){
			
			for(i=0;i<symbol;i++){
				Iaxis[i]=0;
				Qaxis[i]=0;
				Idecode[i]=0;
				Qdecode[i]=0;
			}
			
			/*--------------- 情報ビットの発生 -------------------------------------*/
			//I軸の情報ビットの発生//
			for(i=0;i<symbol/code;i++){	//情報ビットを u3,u5,u6,u7 とする
				Iaxis[(i*code)+2]=rand()%2;
				Iaxis[(i*code)+4]=rand()%2;
				Iaxis[(i*code)+5]=rand()%2;
				Iaxis[(i*code)+6]=rand()%2;
			}
			
			//Q軸の情報ビットの発生//
			for(i=0;i<symbol/code;i++){	//情報ビットを u3,u5,u6,u7 とする
				Qaxis[(i*code)+2]=rand()%2;
				Qaxis[(i*code)+4]=rand()%2;
				Qaxis[(i*code)+5]=rand()%2;
				Qaxis[(i*code)+6]=rand()%2;
			}
			
			/*-----------------------------------------------------------------*/
			
			
			
			
			/*--------------- 符号語の作成 -------------------------------------*/
			//I軸の情報ビットの作成//
			for(i=0;i<symbol/code;i++){
				Iaxis[(i*code)+0]=(Iaxis[(i*code)+2]+Iaxis[(i*code)+4]+Iaxis[(i*code)+6])%2;	//α=u4+u5+u6+u7 であるので u4 は、α が0または1で値が変わる
				Iaxis[(i*code)+1]=(Iaxis[(i*code)+2]+Iaxis[(i*code)+5]+Iaxis[(i*code)+6])%2;
				Iaxis[(i*code)+3]=(Iaxis[(i*code)+4]+Iaxis[(i*code)+5]+Iaxis[(i*code)+6])%2;
			}
			
			//Q軸の情報ビットの作成//
			for(i=0;i<symbol/code;i++){
				Qaxis[(i*code)+0]=(Qaxis[(i*code)+2]+Qaxis[(i*code)+4]+Qaxis[(i*code)+6])%2;
				Qaxis[(i*code)+1]=(Qaxis[(i*code)+2]+Qaxis[(i*code)+5]+Qaxis[(i*code)+6])%2;
				Qaxis[(i*code)+3]=(Qaxis[(i*code)+4]+Qaxis[(i*code)+5]+Qaxis[(i*code)+6])%2;
			}
			
			/*-----------------------------------------------------------------*/
			// for(i=0;i<symbol/code;i++){	//確認 0
				// for(j=0;j<code;j++){
				// printf("I[%d]=%d\t Q[%d]=%d\n",j,Iaxis[j],j,Qaxis[j]);
				// }
			// }
			
			
			
			
			/*--------------- QPSKベースバンド信号の発生 --------------------------*/
			// ベースバンド信号の初期化//
			// for(i=0;i<symbol/code;i++){
				// base[(i*code)+1]=0.0;
			// }
			// 何もなし
			
			for(i=0;i<symbol;i++){	//シンボル数ぶんのベースバンド信号の振幅決定
				if (Iaxis[i]==0){
					if (Qaxis[i]==0)
						base[i]=1/sqrt(2*s)+1/sqrt(2*s)*I;
					else //Qaxis[i]=1のとき
						base[i]=1/sqrt(2*s)-1/sqrt(2*s)*I;
				}
				else {//Iaxis[i]=1のとき
					if (Qaxis[i]==0)
						base[i]=-1/sqrt(2*s)+1/sqrt(2*s)*I;
					else //Qaxis[i]=1のとき
						base[i]=-1/sqrt(2*s)-1/sqrt(2*s)*I;
				}
			}
			
			for(i=0;i<symbol;i++){		//ベースバンド信号をシンボル数からサンプル数に再格納
				// for(j=s*(i-1);j<s*(i);j++){
					// base[j]+=base[i];
				// }
				for(j=s*i;j<s*(i+1);j++){
					basedash[j]=base[i];
				}
			}
			
			/*-----------------------------------------------------------------*/
			// for(i=0;i<symbol;i++){	//確認 0
				// printf("%f+%f\n",creal(base[i]),cimag(base[i]));
			// }
			// for(i=0;i<N;i++){	//確認 0
				// printf("%f+%f\n", creal(basedash[i]),cimag(basedash[i]));
			// }
			
			
			
			
			/*--------------- 複素ガウス雑音の発生 -------------------------------*/
			sigma=sqrt(pow(10.0,-0.1*ENdB)*BT/(2.0*s*BIT)*code/data);
			AWGN(sigma,gauss);
			
			//QPSK信号へのガウス雑音の付加//
			for(i=0;i<N;i++){
				basegauss[i]=basedash[i]+gauss[i];
			}
			
			/*-----------------------------------------------------------------*/
			// for(i=0;i<N;i++){	//確認 0
				// printf("%f+%f\n",creal(gauss[i]),cimag(gauss[i]));
			// }
			// for(i=0;i<N;i++){	//確認 0
				// printf("%f\t %f\n",creal(basegauss[i]),cimag(basegauss[i]));
			// }
			
			
			
			
			/*--------------- 復調 -----------------------------------------------------------------------------------------------------------------------------*/
			for(i=0;i<symbol;i++){
				if(creal(basegauss[i*s])<0){
					Iaxdecode[i]=1;
				}
				else if(creal(basegauss[i*s])>=0){
					Iaxdecode[i]=0;
				}
				if(cimag(basegauss[i*s])<0){
					Qaxdecode[i]=1;
				}
				else if(cimag(basegauss[i*s])>=0){
					Qaxdecode[i]=0;
				}
			}
			
		
			//これいる？
			for(i=0;i<symbol;i++){
				Idecode[i]=Iaxdecode[i];
				Qdecode[i]=Qaxdecode[i];
			}
			
						// for(i=0;i<symbol;i++){	//確認
							// printf("creal(basegauss[%d])=%f\t Idecode[%d]=%d\n",i,creal(basegauss[i]),i,Idecode[i]);
						// }
						// for(i=0;i<symbol;i++){	//確認
							// printf("creal(basegauss[%d])=%f\t Qdecode[%d]=%d\n",i,creal(basegauss[i]),i,Qdecode[i]);
						// }
															
			for(i=0;i<symbol/code;i++){		//for(i=1;i<=symbol;i++)のときIdecode[((i-1)*code)+4]
				//I軸//
				syndrome[0]=(Idecode[(i*code)+3] + Idecode[(i*code)+4] + Idecode[(i*code)+5] + Idecode[(i*code)+6])%2;
				syndrome[1]=(Idecode[(i*code)+1] + Idecode[(i*code)+2] + Idecode[(i*code)+5] + Idecode[(i*code)+6])%2;
				syndrome[2]=(Idecode[(i*code)+0] + Idecode[(i*code)+2] + Idecode[(i*code)+4] + Idecode[(i*code)+6])%2;
				
						// printf("s[0]=(%d+%d+%d+%d)%2=%d\n",Idecode[(i*code)+3],Idecode[(i*code)+4],Idecode[(i*code)+5],Idecode[(i*code)+6],syndrome[0]);
						// printf("s[1]=(%d+%d+%d+%d)%2=%d\n",Idecode[(i*code)+1],Idecode[(i*code)+2],Idecode[(i*code)+5],Idecode[(i*code)+6],syndrome[1]);
						// printf("s[2]=(%d+%d+%d+%d)%2=%d\n",Idecode[(i*code)+0],Idecode[(i*code)+2],Idecode[(i*code)+4],Idecode[(i*code)+6],syndrome[2]);
						// printf("%d\n",syndrome[0]*4 + syndrome[1]*2 + syndrome[2]);
				
				if((syndrome[0]*4 + syndrome[1]*2 + syndrome[2])!=0){	//syndrome[0],[1]は偶数項*4,*2がかかっていることから常に偶数。合計の判定基準はsyndrome[2]のみの判定となる。
																		//よって、if条件文を満たす(＝条件式s[0]*4+s[1]*2+s[0]!=1となる）ときはsyndrome[2]=1のとき
					Idecode[i*code + ((syndrome[0]*4 + syndrome[1]*2 + syndrome[2])-1)]=(Idecode[i*code + ((syndrome[0]*4 + syndrome[1]*2 + syndrome[2])-1)]+1)%2;
					// printf("syndrome[0]*4 + [1]*2 [2]=%d\t",syndrome[0]*4 + syndrome[1]*2 + syndrome[2]);
					// printf("Idecode[]=%d\n",Idecode[i*code + ((syndrome[0]*4 + syndrome[1]*2 + syndrome[2])-1)]);
				}//(syndorome[0] syndrome[1] syndrome[2])⇒(2^2 2^1 2^0)
														  //↑ 0 or 1　の2ビットの１乗
				
				//Q軸//
				syndrome[0]=(Qdecode[(i*code)+3] + Qdecode[(i*code)+4] + Qdecode[(i*code)+5] + Qdecode[(i*code)+6])%2;
				syndrome[1]=(Qdecode[(i*code)+1] + Qdecode[(i*code)+2] + Qdecode[(i*code)+5] + Qdecode[(i*code)+6])%2;
				syndrome[2]=(Qdecode[(i*code)+0] + Qdecode[(i*code)+2] + Qdecode[(i*code)+4] + Qdecode[(i*code)+6])%2;
						// printf("s[0]=(%d+%d+%d+%d)%2=%d\n",Qdecode[(i*code)+3],Qdecode[(i*code)+4],Qdecode[(i*code)+5],Qdecode[(i*code)+6],syndrome[0]);
						// printf("s[1]=(%d+%d+%d+%d)%2=%d\n",Qdecode[(i*code)+1],Qdecode[(i*code)+2],Qdecode[(i*code)+5],Qdecode[(i*code)+6],syndrome[1]);
						// printf("s[2]=(%d+%d+%d+%d)%2=%d\n",Qdecode[(i*code)+0],Qdecode[(i*code)+2],Qdecode[(i*code)+4],Qdecode[(i*code)+6],syndrome[2]);
						// printf("%d\n",syndrome[0]*4 + syndrome[1]*2 + syndrome[2]);

				if((syndrome[0]*4 + syndrome[1]*2 + syndrome[2])!=0){
					Qdecode[i*code + ((syndrome[0]*4 + syndrome[1]*2 + syndrome[2])-1)]=(Qdecode[i*code + ((syndrome[0]*4 + syndrome[1]*2 + syndrome[2])-1)]+1)%2;
					// printf("syndrome[0]*4 + [1]*2 [2]=%d\t",syndrome[0]*4 + syndrome[1]*2 + syndrome[2]);
					// printf("Qdecode[]=%d\n",Qdecode[i*code + ((syndrome[0]*4 + syndrome[1]*2 + syndrome[2])-1)]);
				}
			}
			/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
			
			
			
			
			/*--------------- 誤り判定 -----------------------------------------*/
			for(i=0;i<symbol;i++){
				icount+=(Iaxis[i]+Idecode[i])%2;	//(Iaxis[i],Idecode[i])=(0,0) and (1,1)のときはicount(or qcount)は %2 の結果、0になる。
				qcount+=(Qaxis[i]+Qdecode[i])%2;	//一方、(Iaxis[i],Idecode[i])=(1,0) and (0,1)のときはicount(or qcount)は %2 の結果、1になり、カウントが加算されていく。
			}
			
			if(icount+qcount<1000){
				loop++;
			}
		}
			// BER=(double)(icount+qcount)/(BIT*loop*symbol);
			// printf("%e\n",BER);
			

		BER=(double)(icount+qcount)/(BIT*loop*symbol);
		printf("%d\t", loop);
		printf("%d\t %e\n",ENdB,BER);
	}
	
	
	
	
	return (0);
}



/*================================== 関数部 ===========================================*/
void DFT(double complex Freq[], double complex Time[]){//DFT（離散フーリエ変換）の作成//

	int k,n;
	double shoulder;
	for(k=0;k<N;k++){
		Freq[k]=0;
	}
	for(k=0;k<N;k++){
		for(n=0; n<N; n++){
			shoulder=2*M_PI*(k-(N/2))*(n-(N/2))/N;		//なぜかN/2をk,nから引くらしい・・・。
			Freq[k]+=(1/sqrt(N))*Time[n]*cexp(-shoulder*I);
		}
	}
}

void IDFT(double complex Time[],double complex Freq[]){//IDFT（逆離散フーリエ変換）の作成//

	int k,n;
	double shoulder;
	for(n=0;n<N;n++){
		Time[n]=0;
	}
	for(n=0;n<N;n++){
		for(k=0;k<N;k++){
			shoulder=2*M_PI*(k-(N/2))*(n-(N/2))/N;		//なぜかN/2をk,nから引くらしい・・・。
			Time[n]+=(1/sqrt(N))*Freq[k]*cexp(shoulder*I);
			}
	}
}

void FFT(double complex Time[],double complex Freq[]){
	
	int i,j,k,n;
	double complex x[N],y[N],w[N];

	for(n=0;n<N;n++){
		x[n]=Time[n];
	}
	
	
	
	for(i=2;i<=N;i*=2){
	
		
		for(j=0;j<N/i;j++){
			
			for(k=0;k<i/2;k++){
			
				y[i*j+k]   = x[i*j/2+k] + cexp((-2*M_PI*I*(k-N/2))/i)*x[i*j/2+k+N/2];
				y[i*j+k+i/2] = x[i*j/2+k] - cexp((-2*M_PI*I*(k-N/2))/i)*x[i*j/2+k+N/2];
				
			}
		
		}
		
		for(n=0; n<N; n++){
	
		x[n]=y[n];
		
	}	
		for(n=0;n<N;n++){
		y[n]=0.0+0.0*I;
	}
		
		
	}
	

	
	
	for(n=0; n<N; n++){
	
	Freq[n]=x[n];
	Freq[n]/=sqrt(N);
	
	}

}	

void IFFT(double complex Freq[],double complex Time[]){
	
	int i,j,k,n;
	double complex x[N],y[N];

	for(n=0;n<N/2;n++){
		x[n]=Freq[n+N/2];
		x[n+N/2]=Freq[n];
	}

	for(i=2;i<=N;i*=2){
	
		
		for(j=0;j<N/i;j++){
			
			for(k=0;k<i/2;k++){
			
				y[i*j+k]   = x[i*j/2+k] + cexp((2*M_PI*I*k)/i)*x[i*j/2+k+N/2];
				y[i*j+k+i/2] = x[i*j/2+k] - cexp((2*M_PI*I*k)/i)*x[i*j/2+k+N/2];
				
			}
		
		}
		
		for(n=0; n<N; n++){
	
		x[n]=y[n];
		
	}	
		for(n=0;n<N;n++){
		y[n]=0.0+0.0*I;
	}
		
		
	}
	
	for(n=0; n<N; n++){
	
	Time[n]=x[n];
	Time[n]/=sqrt(N);
	
	}

}	

void QPSKbaseband(double complex base[]){	//QPSKベースバンド信号//

	int n;
	int bi,bq;
	int Bi[symbol],Bq[symbol];
	
	for(n=0;n<symbol;n++){
		bi=rand()%2;	//I軸の0,1の発生
		bq=rand()%2;	//Q軸の0,1の発生
		Bi[n]=bi;		//I軸の数値格納
		Bq[n]=bq;		//Q軸の数値格納
	}
	
	for(n=0;n<symbol;n++){	//シンボル数ぶんのベースバンド信号の振幅決定
		if (Bi[n]==0){
			if (Bq[n]==0)
				base[n]=1/sqrt(2*s)+1/sqrt(2*s)*I;
			else //Bq[n]=1のとき
				base[n]=1/sqrt(2*s)-1/sqrt(2*s)*I;
		}
		else {//Bi[n]=1のとき
			if (Bq[n]==0)
				base[n]=-1/sqrt(2*s)+1/sqrt(2*s)*I;
			else //Bq[n]=1のとき
				base[n]=-1/sqrt(2*s)-1/sqrt(2*s)*I;
		}
	}
}

void SIXTQAMbaseband(double complex base[]){	//16QAM(SiXTeenQAM)ベースバンド信号//

	int n;
	int bi1,bi2,bq1,bq2;
	int Bi1[symbol],Bi2[symbol],Bq1[symbol],Bq2[symbol];
	double baseamp;
		
	for(n=0;n<symbol;n++){	//I軸の2ビット分の発生//
		bi1=rand()%2;
		bi2=rand()%2;
		Bi1[n]=bi1;
		Bi2[n]=bi2;
	}
	for(n=0;n<symbol;n++){	//Q軸の2ビット分の発生//
		bq1=rand()%2;
		bq2=rand()%2;
		Bq1[n]=bq1;
		Bq2[n]=bq2;
	}
	
	baseamp=1/sqrt(10*s);		//ベースバンド信号の基本振幅
	
	for(n=0;n<symbol;n++){
		base[n]=0;		//Baseの初期化
			if (Bi1[n]==0 && Bi2[n]==0){			//base[n]=a+jbで表される。まずはa(=I軸)のみを求める
				base[n]+=3*baseamp;
			}
			else if (Bi1[n]==0 && Bi2[n]==1){
				base[n]+=baseamp;
			}
			else if (Bi1[n]==1 && Bi2[n]==1){
				base[n]+=-baseamp;
			}
			else if (Bi1[n]==1 && Bi2[n]==0){
				base[n]+=-3*baseamp;
			}										
			
			if (Bq1[n]==0 && Bq2[n]==0){
				base[n]+=3*baseamp*I;
			}
			else if (Bq1[n]==0 && Bq2[n]==1){
				base[n]+=baseamp*I;
			}
			else if (Bq1[n]==1 && Bq2[n]==1){
				base[n]+=-baseamp*I;
			}
			else if (Bq1[n]==1 && Bq2[n]==0){
				base[n]+=-3*baseamp*I;
			}					
	}
}

void AWGN(double sigma,double complex cgn[]){	//AWGNの作成//
	int i;
	double Ix,Qy,z1,z2;
		
	for(i=0;i<N;i++){
		Ix=((double)rand()+1.0)/((double)RAND_MAX+1.0);	//一様乱数
		Qy=((double)rand()+1.0)/((double)RAND_MAX+1.0);	//一様乱数
		z1=sqrt(-2.0*log(Ix))*cos(2.0*M_PI*Qy)*(sigma);	//ガウス分布におけるI軸の時間波形
		z2=sqrt(-2.0*log(Ix))*sin(2.0*M_PI*Qy)*(sigma);	//ガウス分布におけるQ軸の時間波形
		cgn[i]=z1+z2*I;				//複素ガウス雑音波形 Complex Gaussian Noise//
		//printf("%f\t%f\n",Ix,Qy);
	}
}