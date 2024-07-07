////////////  OFDM-QPSK,OFDM-16QPSKのPAPR特性の算出 及び SC-QPSK,SC-16QAMのPAPR特性の算出 /////////
/*
・SC伝送はルートロールオフ伝送(ルートロールオフ率α=0.5)におけるPAPR特性

*/
#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム

#define a 1						//波数
#define s 8						//オーバーサンプル比
#define K 256					//サブキャリア数
#define N (s*K)					//全サンプル数
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1シンボル時間
#define fc a/Ts 				//搬送波周波数fc
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

///　root-Nyquist LPF
#define alpha 0.5				//ロールオフ率（Roll-off Ratio）
#define fs 1/Ts					//サンプリング周波数

/// PAPR特性，確率密度関数
#define loop 100000					//試行回数
#define range 1500				//試行範囲
#define w 0.01					//範囲の幅

////☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆目的☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆☆////
// (1)QPSK-OFDM，(2)16QAM-OFDM (3)QPSK-SC (4)16QAM-SC　のPAPR特性算出する ///////


////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void QPSK(double complex base[]);
void QAM16(double complex QAM[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void root_Nyquist(double root_ratio,double complex Gt[],double complex Gr[]);
void SP(double complex Serial[],double complex Parallel[]);
void PS(double complex Parallel[],double complex Serial[]);
double max2(double complex data[]);	//データの2乗の最大値を返す関数
double average2(double complex data[]);	//データの2乗(電力)の平均値を返す

int main(void)
{	
	int i,j,m,k; 
	
///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 

////////////////////// 時間と周波数の設定　　//////////////////////
	double t[N]; 						//時間t    0～K*Tsまでの時間
	double f[N];						//周波数f  -Ts/s～Ts/sまでの周波数
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}

///////////////////////////////////　root-Nyquist LPF　構築 ////////////////////////////////////////////////
	double complex Gt[N],Gr[N];		//送受信フィルタの伝達関数

	root_Nyquist(alpha,Gt,Gr);
	
///////////////////////////////////////////////////////////////////////////////////////////////////////
	//// OFDM伝送時に用いる信号 ////////////////////////
	double complex sF1[N];		//一次変調信号
	double complex sF2[N];		//先頭サンプルのみ一次変調信号
	double complex sF[N];		//S/P変換後の信号
	///// シングルキャリア伝送時の信号 /////////////////////
	double complex st1[N];		//QPSK or 16QAM信号
	double complex st2[N];		//先頭サンプルのみ信号
	double complex fft_s2[N];		//st2のフーリエ変換
	double complex sF3[N];		//送信フィルタ通過後の周波数領域信号
	
	/// QPSK,16QAM-OFDMまたはQPSK,16QAM-SC　共通の送信信号　→PAPRを計算するときにstを用いたいから共通にした
	double complex st[N];		//IFFT後の信号	（OFDM送信信号）
	
	///  PAPR　特性算出時のパラメータ /////////////////////////////
	int c1[range];				//OFDMのPAPR用カウンタ
	double PAPR;				//OFDMのPAPR	
	double pdf_PAPR[range];		//OFDMのPAPRの確率密度関数
	double CCDF[range];			//CCDF
	
	for(i=0; i<range; i++)
		c1[i]=0;		//初期化
		
	printf("K=%d \n",K);
	for(m=0;m<loop;m++){	//1回のloopに対して1つのPAPRが求まる
		///// (1)QPSK-OFDM生成 ///////////////////////////////////////////////////
		// QPSK(sF1);			//QPSK信号生成
		// inter(sF1,sF2);		//K本のサブキャリア=先頭サンプルのみ抜き出す
		// SP(sF2,sF);
		// IFFT(sF,st);
		
		///// (2)16QAM-OFDM生成 ///////////////////////////////////////////////////
		// QAM16(sF1);			//QPSK信号生成
		// inter(sF1,sF2);		//K本のサブキャリア=先頭サンプルのみ抜き出す
		// SP(sF2,sF);			
		// IFFT(sF,st);
		
		////// (3)QPSK-SC生成 (4)16QAM-SC生成/////////////////////////////////////////////////////
		
		
		QPSK(st1);			//QPSK信号生成
		// QAM16(st1);		//16QAN信号生成
		inter(st1,st2);		//先頭サンプルのみQPSK信号
		FFT(st2,fft_s2);		//フィルタ通過準備		
		
		for(i=0;i<N;i++)
			sF3[i]=fft_s2[i]*Gt[i];	//送信フィルタ通過
		
		IFFT(sF3,st);				//送信信号完成
		
		///// QPSK_OFDMまたは16QAM_OFDMのPAPR[dB]を算出  //////////////////////////////////////
		PAPR=10*log10(max2(st)/average2(st));
		
		// printf("%f   %fdB\n",max2(st)/average2(st),PAPR);
		
		////////////　　QPSK-OFDMまたは16QAM-OFDM送信信号のPAPRの確率密度関数   //////////////
		for(k=0; k<range; k++){			//PAPRの分類[-0.1,0)，[0,0.1),.....
			if((k-1)*w<=PAPR && PAPR<k*w) 
				c1[k]++;
		}
	}
	
	//　確率密度関数の算出  //////////////////////////////////////////////
	for(k=0; k<range; k++)
		pdf_PAPR[k]=c1[k]/(loop*w);		//データはloop個だから K*loop個ではない
	
	////// CCDF 算出　///////////////////////////////////////////////
	
	
	for(k=0; k<range; k++){		//CCDF=PAPRがある値kを超える確率
		CCDF[k]=0.0;		//初期化
		for(j=k;j<range;j++)	//k以上の確率密度関数を足し合わせる
		CCDF[k]+=w*pdf_PAPR[j];		//幅w×確率密度関数(=確率)の累積
	}
	
	putchar('\n');
	printf("PAPR\t\tPDF\t\tCCDF\t\tK=%d\n",K);
	for(k=0;k<range;k++)
		printf("%f\t%e\t%e\n",k*w,pdf_PAPR[k],CCDF[k]);
	
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

//////////////////////////　ルートオフ率αのルートロールオフ伝送系の構築　/////////////////////////////////////////
void root_Nyquist(double root_ratio,double complex Gt[],double complex Gr[]){

	int i,j;
	double complex *h1;				//ルートロールオフフィルタのインパルス応答
	double *time1,arg1;
	double complex *fft_h1;			//ルートロールオフフィルタのh1(t)の周波数成分	
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
		else if(time1[i]==1/(2*root_ratio*fs) || time1[i]==-1/(2*root_ratio*fs))	
			h1[i]=sqrt((double)N)*root_ratio/2*sin(M_PI/(2*root_ratio));		// 0/0となる部分の極限はπ/4だから
		else
			h1[i]=sqrt((double)N)*sin(arg1)/arg1*cos(root_ratio*arg1)/(1-pow(2*root_ratio*fs*time1[i],2));
	}
//////////////////////////////////　　ロールオフフィルタの周波数応答　　//////////////////////////////////////////
	FFT(h1,fft_h1);
///////////////////////////////////////////　ルート配分　////////////////////////////////////////////////
	for(i=0;i<N;i++){
		Gt[i]=sqrt(cabs(fft_h1[i]))+I*0.0;	//sqrt(double)なので実部，虚部分けてやる
		Gr[i]=Gt[i];	
	}
	
	// putchar('\n');
	// printf("ナイキストロールオフフィルタのh(t)\tH(f)\n");
	// for(i=0;i<N;i++)
		// printf("%f\t\t\t%f\n",cabs(h1[i]),cabs(fft_h1[i]));
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

////////////////////////// データの2乗の最大値を返す関数 /////////////////////////
double max2(double complex data[]){
	int i;
	double max;
	
	max=cabs(data[0])*cabs(data[0]);
	for(i=1;i<N;i++){
		if(cabs(data[i])*cabs(data[i])>max)
			max=cabs(data[i])*cabs(data[i]);
	}
	
	return max;
}
	
double average2(double complex data[]){
	int i;
	double sum=0.0;
	
	for(i=0;i<N;i++)
		sum+=cabs(data[i])*cabs(data[i]);
	
	return sum/N;
}
	










