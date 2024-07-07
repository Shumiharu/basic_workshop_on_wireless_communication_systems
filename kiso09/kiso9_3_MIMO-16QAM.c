/*/////////////////////  多重伝送 MIMO-OFDM-16QAM　ZF，MMSE，MLD等化　　///////////////////////////////
・MLD_16QAMの配列がQ^Nt乗でNt=4にしたらコアダンプを起こす→limit -s unlimited　で動的確保を行う
・雑音電力の設定係数(Nt*Nt)/Nrが正しい！

*/
#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム


#define a 1			 			//波数
#define s 1						//オーバーサンプル比
#define K 64					//サブキャリア数
#define N (s*K)					//全サンプル数
#define sigmaN 0.5				//ガウス雑音の分散
#define Ts 0.001				//1サンプル時間(OFDMでは)
#define fc a/Ts 				//搬送波周波数fc
#define T_sample (Ts/s) 		//1サンプル当たりの時間
#define delta_f (1/(Ts*K))		//Δf=1/T Δf:スペクトルの幅  T(=Ts*K):時間長

//OFDM特有のパラメータ↓↓↓↓
#define L 1						//L波マルチパス
#define dB 1					//何dB減衰
#define GI 4					//ガードインターバルのサンプル数

//MIMO特有のパラメータ
#define Nt 2					//送信アンテナNt
#define Nr 2					//受信アンテナNr
#define Q 16 				//変調多値数
#define choise 1			//1=ZF等化，2=MMSE等化，3=MLD

//BER特性用↓↓↓↓
#define bit 4.0 				//1シンボル当たりのビット数 16QAMの場合は4
#define BnTs 1   				
#define first 50000				//いったん雑にBERを算出する


////////////////////////関数プロトタイプ宣言/////////////////////////////////
void FFT(double complex Time[],double complex Freq[]);
void IFFT(double complex Freq[],double complex Time[]);
void QPSK(double complex base[]);
void QAM16(double complex QAM[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void ht_MPF(double complex ht_mpf[Nr][Nt][N]);
void SP(double complex Serial[],double complex Parallel[]);
void PS(double complex Parallel[],double complex Serial[]);
void GI_insertion(double complex data[],double complex data_GI[]);
void GI_removal(double complex data_GI[],double complex data[]);
void matrix_prod(int n1,int n2,int n3,int nd,double complex mx1[n1][n2][nd],double complex mx2[n2][n3][nd],double complex mxout[n1][n3][nd]);
void inv_matrix(int num,int nd,double complex IN[num][num][nd],double complex OUT[num][num][nd]);
void Hermitian_transpose(int n1,int n2,int nd,double complex mxin[n1][n2][nd],double complex mxout[n2][n1][nd]);
void ZF(double complex H[Nr][Nt][N],double complex r1[Nr][N],double complex r2[Nt][N]);
void MMSE(double complex H[Nr][Nt][N],double variance,double complex r1[Nr][N],double complex r2[Nt][N]);
void MLD_QAM16(double complex H[Nr][Nt][N],double complex rt[Nr][N],double complex estimate[Nt][N]);

//////0～n-1まで数を重複を許して並べる順列をpermに格納///////
void perm_sub(int perm[],int n1,int n2,int n3[]);
void copy_perm(int perm1[],int num[]);

int main(void)
{	
	int i,j,k,n,m,l; 
	
///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 

////////////////////// 時間と周波数の設定　　//////////////////////
	double t[N]; 						//時間t    0～K*Tsまでの時間
	double f[N];						//周波数f  -Ts/s～Ts/sまでの周波数
	
	for(i=0;i<N;i++){
		t[i]=T_sample*i;
		f[i]=delta_f*(i-N/2);
	}

////////////////////////
	
//////////////////↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓送信側↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓//////////////////////////////	
	double amp=1/sqrt(10.0*s);
	double complex (*sF1)[N]=(double complex(*)[N])malloc(Nt*N*sizeof(double complex));		//QPSK信号
	double complex (*sF2)[N]=(double complex(*)[N])malloc(Nt*N*sizeof(double complex));		//先頭サンプルのみQPSK信号
	double complex (*sF)[N]=(double complex(*)[N])malloc(Nt*N*sizeof(double complex));		//S/P変換後の信号
	double complex (*st)[N]=(double complex(*)[N])malloc(Nt*N*sizeof(double complex));		//IFFT後時間信号（雑音状）
	double complex (*st_GI)[(K+GI)*s]=(double complex(*)[(K+GI)*s])malloc(Nt*((K+GI)*s)*sizeof(double complex));	//GI挿入後の時間信号
	double complex (*sF_GI)[(K+GI)*s]=(double complex(*)[(K+GI)*s])malloc(Nt*((K+GI)*s)*sizeof(double complex));	//GI挿入後の周波数領域信号
	double complex (*tmp)[N]=(double complex(*)[N])malloc(Nt*N*sizeof(double complex));		//次サブキャリアへの干渉のため自シンボルのデータを格納
	
///////////////////↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓通信路行列生成用↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓//////////////////////////////
	double complex (*ht_multifd)[Nt][N]=(double complex(*)[Nt][N])malloc(Nr*Nt*N*sizeof(double complex));		//マルチパスフェージングのインパルス応答
	double complex (*Hf_multifd)[Nt][N]=(double complex(*)[Nt][N])malloc(Nr*Nt*N*sizeof(double complex));		//マルチパスフェージングの伝達関数
	
	
//////////////////↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓受信側↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓//////////////////////////////
	double complex (*st3)[(K+GI)*s]=(double complex(*)[(K+GI)*s])malloc(Nr*((K+GI)*s)*sizeof(double complex));	//マルチパス通信路通過後の各受信アンテナでの受信信号
	double complex (*Delay)[(K+GI)*s]=(double complex(*)[(K+GI)*s])malloc(Nr*((K+GI)*s)*sizeof(double complex));	//受信アンテナ0,1,...Nrに到来する遅延波
	double complex (*awgn)[(K+GI)*s]=(double complex(*)[(K+GI)*s])malloc(Nr*((K+GI)*s)*sizeof(double complex));		//ガウス雑音
	double complex (*rt1)[(K+GI)*s]=(double complex(*)[(K+GI)*s])malloc(Nr*((K+GI)*s)*sizeof(double complex));		//ガウス雑音重畳後の時間信号
	double complex (*rt2)[N]=(double complex(*)[N])malloc(Nr*N*sizeof(double complex));		//GI除去後の時間信号
	double complex (*rF2)[N]=(double complex(*)[N])malloc(Nr*N*sizeof(double complex));		//rt2のフーリエ変換
	
//////////////////↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓受信側判定時(要素数がNtになる)↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓//////////////////////////////
	double complex (*rhatF)[N]=(double complex(*)[N])malloc(Nt*N*sizeof(double complex));	//FDE後の周波数信号
	double complex (*rF)[N]=(double complex(*)[N])malloc(Nt*N*sizeof(double complex));		//rhatをP/S変換

///////////////////////////////////BER特性算出時に用いる変数　　/////////////////////////////////////////////////
	double sigma;							//ガウス雑音の分散			
	int EbN0;								//Eb/N0
	double ci,cq;							//BERカウント用
	double ber;  							//誤り率
	int loop;
	int LOOP;			//雑にいったんBERを出す用の試行回数
	
	
	//等化方法を表示
	if(choise==1)	printf("ZF等化\n");
	if(choise==2)	printf("MMSE等化\n");	
	if(choise==3)	printf("MLD\n");	
	
	
	printf("%d×%dMIMO-OFDM(16QAM)%dパス準静的レイリーフェージング，%ddB減衰モデル　GI=%d\n",Nt,Nr,L,dB,GI);
	printf("ρ[dB]        BER \n");
	for(EbN0=0; EbN0<=50; EbN0+=5){
		ci=0.0; cq=0.0; 
		LOOP=first;			//Eb/N0の値ごとにfirst回でBERを軽く算出する
		
		//↓↓↓↓↓↓↓雑音電力の設定注意↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
		sigma=(Nt*Nt)/Nr*BnTs/(2*s*bit)*pow(10.0,(-0.1*EbN0))*((K+GI)/(double)K);		//GIの挿入によって帯域幅がK+GI/K倍拡大→雑音電力K+GI/K倍される
		//  遅延波格納用tmp　初期化 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////　　
		 
		for(i=0;i<N;i++){
			for(j=0;j<Nt;j++)
			tmp[j][i]=0+0*I;
		}
		
		/// 　フレーム数分OFDMを伝送する　//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for(loop=0; loop<LOOP; loop++){	
			//// 各送信アンテナでフレームごとに送信信号stを作成 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for(j=0;j<Nt;j++){	//送信アンテナの数だけデータを生成
				QAM16(sF1[j]);						//QAM16信号生成
				inter(sF1[j],sF2[j]);					//K本のサブキャリア=先頭サンプルのみ抜き出す
				SP(sF2[j],sF[j]);						//S/P変換，周波数軸の中央に配置
				IFFT(sF[j],st[j]);		
				GI_insertion(st[j],st_GI[j]);			//GI挿入
			}
			
			ht_MPF(ht_multifd);			//マルチパス通信路のインパルス応答		
			
			for(i=0;i<Nr;i++){
				for(j=0;j<Nt;j++)
					FFT(ht_multifd[i][j],Hf_multifd[i][j]);		//マルチパスフェージングの伝達関数H(f)
			}
		
		
			// 正規化    　※総電力で正規化したら間違った振幅スペクトルが得られてしまった
			for(i=0;i<Nr;i++){
				for(j=0;j<Nt;j++){
					for(k=0;k<N;k++)
						Hf_multifd[i][j][k]=Hf_multifd[i][j][k]/sqrt(1.0/(double)N);		//はい？	
				}
			}				
				
			/// st_GIがマルチパス通信路を通過　///////////////////////////////////////////////////////////////
			for(j=0;j<Nr;j++){		//☆☆☆
				for(i=0; i<(K+GI)*s; i++)			
					st3[j][i]=0.0+0.0*I;		//初期化
			}
			
			//畳み込み和（正確）
			for(m=0;m<Nr;m++){			//受信アンテナ0,受信アンテナ1...の順に
				for(k=0;k<Nt;k++){		//送信アンテナ0,送信アンテナ1...の順に
					for(i=0; i<L; i++){		//直接波1個＋遅延波(L-1)個
						for(j=0; j<(K+GI)*s; j++){		
							if(j<i*s)	//i=0の時は直接波なのでこのifには入らない，遅延波1なら1サブキャリア分前回の送信信号tmpをs点コピー
								Delay[m][j]=tmp[k][j+N-(i*s)];	//受信アンテナ0,受信アンテナ1...に到来する遅延波の順に	
							else							
								Delay[m][j]=st_GI[k][j-(i*s)];	//前シンボルの遅延をコピーし残りをst_GIをst_GI[0]から順番にコピーしていく		
								
							st3[m][j]+=Delay[m][j]*ht_multifd[m][k][i*s];	//送信アンテナ0から到来したDelay,送信アンテナ1から到来したDelayを受信アンテナ0で畳み込み和をとる
						
						}			//Delayは送信アンテナ番号kによって中身は違う
					}//☆☆☆Delay
				}
			}

			
			//　次シンボル遅延波用に送信信号をコピー　//////////////////////////////////////////////////////////////
			for(k=0;k<Nt;k++){
				for(i=0; i<N; i++)
					tmp[k][i]=st[k][i];		//次シンボルの遅延波に含まれるのでtmpに一時的に送信信号を格納
			}
			
			///////　ガウス雑音付加　　////////////////////////////////////////////////////////////////////////
			for(m=0;m<Nr;m++)
				boxmuller((K+GI)*s,awgn[m],0.0,sigma);		//　AWGNの生成　
			
			///////　それぞれのアンテナに異なるガウス雑音を加える ///////////////////////////////////////////////////////
			for(m=0;m<Nr;m++){
				for(i=0;i<(K+GI)*s;i++){
					rt1[m][i]=st3[m][i]+awgn[m][i];		//rt(t)=st(t)+AWGN(t)　
					// rt1[m][i]=st3[m][i]+0;			//雑音を入れないver→BER=0となっているか確認用
				}
			}
			
			
			////// GI removal　 ///////////////////////////////////////
			for(m=0;m<Nr;m++)
				GI_removal(rt1[m],rt2[m]);			
			//☆☆☆
			
			////// FFT r2(t)→R2(f)　/////////////////////////////////////////////////////////////
			for(m=0;m<Nr;m++)
				FFT(rt2[m],rF2[m]);
			
			////// FDE(ZF criterion) Rhat(f)=R2(f)/H(f) /////////////////////////////////////////////////////////////
			//等化方法を選択
			if(choise==1)	ZF(Hf_multifd,rF2,rhatF);
			if(choise==2)	MMSE(Hf_multifd,sigma,rF2,rhatF);
			if(choise==3)	MLD_QAM16(Hf_multifd,rF2,rhatF);
	
			//////  P/S変換　Rhat(f)→R(f)　///////////////////////////////////////////////////////////////////////////  
			for(m=0;m<Nt;m++)
				PS(rhatF[m],rF[m]);		//☆☆
		
			//// 誤り判定(先頭点のみ)  ////////////////////////////////////////////////////////////////////////
			for(m=0;m<Nt;m++){	//☆☆☆
				
				
				for(j=0;j<K;j++){
					///////////////////////////  I軸判定　　//////////////////////////////////////////
					////////////    (00)伝送時の判定   /////////////
					if(creal(sF2[m][j*s])>=2*amp){
						if(0<=creal(rF[m][j*s]) && creal(rF[m][j*s])<2*amp){	
							ci++;	//(00)伝送時(01)と誤る
						}else if(-2*amp<=creal(rF[m][j*s]) && creal(rF[m][j*s])<0){
							ci+=2;	//(00)伝送時(11)と誤る
						}else if(creal(rF[m][j*s])<=-2*amp){
							ci++;	//(00)伝送時(10)と誤る
						}
					}
						
					////////////    (01)伝送時の判定   /////////////
					if(0<=creal(sF2[m][j*s]) && creal(sF2[m][j*s])<2*amp){
						if(2*amp<=creal(rF[m][j*s])){	
							ci++;	//(01)伝送時(00)と誤る
						}else if(-2*amp<=creal(rF[m][j*s]) && creal(rF[m][j*s])<0){
							ci++;	//(01)伝送時(11)と誤る
						}else if(creal(rF[m][j*s])<=-2*amp){
							ci+=2;	//(01)伝送時(10)と誤る
						}	
					}
						
					////////////    (11)伝送時の判定   /////////////
					if(-2*amp<=creal(sF2[m][j*s]) && creal(sF2[m][j*s])<0){
						if(2*amp<=creal(rF[m][j*s])){	
							ci+=2;	//(11)伝送時(00)と誤る
						}else if(0<=creal(rF[m][j*s]) && creal(rF[m][j*s])<2*amp){
							ci++;	//(11)伝送時(01)と誤る
						}else if(creal(rF[m][j*s])<=-2*amp){
							ci++;	//(11)伝送時(10)と誤る
						}	
					}
						
					////////////    (10)伝送時の判定   /////////////
					if(creal(sF2[m][j*s])<=-2*amp){
						if(2*amp<=creal(rF[m][j*s])){	
							ci++;	//(10)伝送時(00)と誤る
						}else if(0<=creal(rF[m][j*s]) && creal(rF[m][j*s])<2*amp){
							ci+=2;	//(10)伝送時(01)と誤る
						}else if(-2*amp<=creal(rF[m][j*s]) && creal(rF[m][j*s])<0){
							ci++;	//(10)伝送時(11)と誤る
						}	
					}
						
					///////////////////////////  Q軸判定　　//////////////////////////////////////////
					////////////    (00)伝送時の判定   /////////////
					if(cimag(sF2[m][j*s])>=2*amp){
						if(0<=cimag(rF[m][j*s]) && cimag(rF[m][j*s])<2*amp){	
							cq++;	//(00)伝送時(01)と誤る
						}else if(-2*amp<=cimag(rF[m][j*s]) && cimag(rF[m][j*s])<0){
							cq+=2;	//(00)伝送時(11)と誤る
						}else if(cimag(rF[m][j*s])<=-2*amp){
							cq++;	//(00)伝送時(10)と誤る
						}
					}
						
					////////////    (01)伝送時の判定   /////////////
					if(0<=cimag(sF2[m][j*s]) && cimag(sF2[m][j*s])<2*amp){
						if(2*amp<=cimag(rF[m][j*s])){	
							cq++;	//(01)伝送時(00)と誤る
						}else if(-2*amp<=cimag(rF[m][j*s]) && cimag(rF[m][j*s])<0){
							cq++;	//(01)伝送時(11)と誤る
						}else if(cimag(rF[m][j*s])<=-2*amp){
							cq+=2;	//(01)伝送時(10)と誤る
						}	
					}
						
					////////////    (11)伝送時の判定   /////////////
					if(-2*amp<=cimag(sF2[m][j*s]) && cimag(sF2[m][j*s])<0){
						if(2*amp<=cimag(rF[m][j*s])){	
							cq+=2;	//(11)伝送時(00)と誤る
						}else if(0<=cimag(rF[m][j*s]) && cimag(rF[m][j*s])<2*amp){
							cq++;	//(11)伝送時(01)と誤る
						}else if(cimag(rF[m][j*s])<=-2*amp){
							cq++;	//(11)伝送時(10)と誤る
						}	
					}
						
					////////////    (10)伝送時の判定   /////////////
					if(cimag(sF2[m][j*s])<=-2*amp){
						if(2*amp<=cimag(rF[m][j*s])){	
							cq++;	//(10)伝送時(00)と誤る
						}else if(0<=cimag(rF[m][j*s]) && cimag(rF[m][j*s])<2*amp){
							cq+=2;	//(10)伝送時(01)と誤る
						}else if(-2*amp<=cimag(rF[m][j*s]) && cimag(rF[m][j*s])<0){
							cq++;	//(10)伝送時(11)と誤る
						}	
					}
				}	
			}//誤り判定 end
			
			//BERの値によってループ調整
			if(loop==LOOP-1){
				ber=(ci+cq)/(bit*(double)loop*K*Nt);		//☆☆☆全部でNt個，一旦BERを算出
				
				if(ber>=1E-1 && ber<1)
					LOOP=pow(10,4);			//loop=first=5000→このifに入る→LOOP=10が代入される→最初のforの部分でloop++される→loop<LOOPの判定でloopの方が大きいのでこのforを抜け出せる→BERを計算して結果表示
				if(ber>=1E-2 && ber<1E-1)
					LOOP=pow(10,4);
				if(ber>=1E-3 && ber<1E-2)
					LOOP=pow(10,4);
				if(ber>=1E-4 && ber<1E-3)
					LOOP=pow(10,5);
				if(ber>=1E-5 && ber<1E-4)
					LOOP=5*pow(10,5);
				if(ber>=1E-6 && ber<1E-5)
					LOOP=5*pow(10,6);
				if(ber>=0 && ber<1E-6)
					LOOP=pow(10,7);
				// if(ber>=1E-7 && ber<1E-6)
					// LOOP=pow(10,7);
				// if(ber>=0 && ber<1E-7)
					// LOOP=pow(10,8);
				// printf("BER=%e, loop=%d, LOOP=%d, Eb/N0=%d dB\n",ber,loop,LOOP,EbN0);
			}
			
			ber=(ci+cq)/(bit*(double)(loop+1)*K*Nt);		//☆☆☆☆☆☆
			printf("loop=%d         Eb/N0=%d         %e\r",loop,EbN0,ber);	
			
		}//loop end
		ber=(ci+cq)/(bit*(double)loop*K*Nt);		//☆☆☆全部でNt個
		printf("loop=%d         Eb/N0=%d         %e\n",loop,EbN0,ber);	
	
	}//Eb/N0 end
	
		
			
	
	//メモリ開放
	free(sF1);			free(sF2);				free(sF);
	free(st);			free(st_GI);			free(sF_GI);
	free(tmp);			free(ht_multifd);		free(Hf_multifd);
	free(st3);			free(Delay);			free(awgn);
	free(rt1);			free(rt2);				free(rF2);
	free(rhatF);		free(rF);
	
	
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

///マルチパスフェージング作成(H[0][0][N]={h1,h2,...}，H[0][1][N]={h3,h4,...}，H[1][0][N]={h5,h6,...}，H[1][1][N]={h7,h8,...}) ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ht_MPF(double complex ht_mpf[Nr][Nt][N]){
	int i,j,k,n;
	double P_multifd[L];					//各パスの電力
	double complex awgn1[1];			
	
	P_multifd[0]=(1-pow(10,-(dB/10.0)))/(1-pow(10,-(dB*L/10.0)));
	for(k=0;k<Nr;k++){
		for(n=0;n<Nt;n++){
			for(i=0;i<K;i++){
				for(j=0;j<s;j++){				
					if(i<L && j==0){		//L波分作成　先頭にのみインパルスの値を入れる
						P_multifd[i]=P_multifd[0]*pow(10,-(dB*i/10.0));		//電力を減衰させてゆく P_multifd[0]=P_multifd[0]*1となる．
						boxmuller(1,awgn1,0.0,P_multifd[i]/2.0);			//2σ^2=P_multifd[i]となるようにガウス雑音を作成	
						ht_mpf[k][n][j+s*i]=awgn1[0];
					}else
						ht_mpf[k][n][j+s*i]=0+0*I;
				}
			}
		}
	}//☆☆☆
	
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

///////////////////////GI insertion //////////////////////////////////////
void GI_insertion(double complex data[],double complex data_GI[]){
	int i,j;

	for(j=0;j<(K+GI);j++){	//サブキャリア数+GI長だけ繰り返す
		if(j<GI){		//GI長の数だけこのifの中に入る
			for(i=0;i<s;i++)
				data_GI[i+j*s]=data[s*(K-GI)+j*s+i];	//GI長×s点のデータコピー
		}else{
			for(i=0;i<s;i++)
				data_GI[i+j*s]=data[i+(j-GI)*s];		//data_GIの残りにｓｔ[N]のデータをコピー m=GIとなってるから(m-GI)とする
		}		
	}	
	
	//GI挿入正しいか確認
	/*  putchar('\n');
	printf("番号\tGI挿入後の信号data_GI\t\tGI挿入前の信号data\n");
	for(i=0;i<(K+GI)*s;i++){
		printf("%d番目\t%f\t%f\t",i,creal(data_GI[i]),cimag(data_GI[i]));
		if(i>=(GI*s)){
			printf("%f\t%f\n",creal(data[i-(GI*s)]),cimag(data[i-(GI*s)]));	
		}else{
			printf("%f\t%f\n",0.0,0.0);		//dataの足りない部分を0入れておく
		}
	}  *///これによってdata_GIの最初のGi*s個分dataの後ろのデータがコピーできていることが確認できた．
	
}

////////////////////////// GI removal　 ///////////////////////////////////////
void GI_removal(double complex data_GI[],double complex data[]){
	int i,j;
	
	for(i=0;i<K;i++){
		for(j=0;j<s;j++)
			data[j+i*s]=data_GI[j+(i+GI)*s];		//rt1の最初のGI*s分を捨てて，それ以降のN個のデータをコピーする
	}
			
}

//////////////////////////////////// num×numの行列の逆行列を求める関数 ////////////////////////////

void inv_matrix(int num,int nd,double complex IN[num][num][nd],double complex OUT[num][num][nd]){
	int i,j,k,l;
	double complex temp1,temp2,temp[num][num][nd];
	
	//tempにINを格納
	for(l=0;l<nd;l++){
		for(i=0;i<num;i++){
			for(j=0;j<num;j++){
				temp[i][j][l]=IN[i][j][l];
			}
		}
	}
	
	//単位行列生成
	for(l=0;l<nd;l++){
		for(i=0;i<num;i++){
			for(j=0;j<num;j++){
				if(i==j){
					OUT[i][j][l]=1.0;		
				}
				else{
					OUT[i][j][l]=0.0;
				}
			}
		}
	}
	
	//掃き出し
	for(l=0;l<nd;l++){
		for(i=0;i<num;i++){
			temp1=temp[i][i][l];//i行目でほかの行を割る
			for(j=0;j<num;j++){
				temp[i][j][l]=temp[i][j][l]/temp1;//単位行列を目指す
				OUT[i][j][l]=OUT[i][j][l]/temp1;//逆行列を目指す
			}
			for(k=0;k<num;k++){
				if(i!=k){//左側の行列が単位行列でないとき
					temp2=temp[k][i][l];//ななめの部分の数をひいて0にする
					for(j=0;j<num;j++){
						temp[k][j][l]=temp[k][j][l]-temp[i][j][l]*temp2;
						OUT[k][j][l]=OUT[k][j][l]-OUT[i][j][l]*temp2;
					}
				}
			}
		}
	}
}

/////////////////  num1×num2とnum2×num3行列の掛け算 出力はnum1×num3の行列 ///////////////////////////////// 
//									   	データ数
void matrix_prod(int n1,int n2,int n3,int nd,double complex mx1[n1][n2][nd],double complex mx2[n2][n3][nd],double complex mxout[n1][n3][nd]){
	int i,j,k,l;
	
	for(i=0;i<nd;i++){
		for(j=0;j<n1;j++){
			for(k=0;k<n3;k++){
				mxout[j][k][i]=0.0+0.0*I;
			}
		}
	}
	
	
	
	for(i=0;i<nd;i++){
		for(j=0;j<n1;j++){
			for(k=0;k<n3;k++){
				for(l=0;l<n2;l++){
					mxout[j][k][i]+=mx1[j][l][i]*mx2[l][k][i];
				}
			}
		}
	}
}

//////////////////////////////// エルミート転置行列(Hermitian tranpose)　/////////////////////////////////
// 			　　　　　　		通信路行列　　					等化前信号　　　				 等化後信号　　　
void Hermitian_transpose(int n1,int n2,int nd,double complex mxin[n1][n2][nd],double complex mxout[n2][n1][nd]){
	int i,j,k;
	
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			for(k=0;k<nd;k++)
				mxout[j][i][k]=conj(mxin[i][j][k]);
		}
	}
	
	// エルミート転置になっている確認　両者が符合違いの関係なら正しい
	/* putchar('\n');
	for(i=0;i<nd;i++){
		for(k=0;k<num;k++){
			for(j=0;j<num;j++){
				printf("H=%f,H^H=%f\n",cimag(mxin[k][j][i]),cimag(mxout[j][k][i]));
			}
		}
	} */
	
}

//////////////////////////////// ZF基準　/////////////////////////////////
// 			　　　　　　		通信路行列　　					等化前信号　　　				 等化後信号　　　
void ZF(double complex H[Nr][Nt][N],double complex r1[Nr][N],double complex r2[Nt][N]){
	int i,j,k;
	double complex H_inv[Nr][Nt][N];
	double complex tmp1[Nr][1][N],tmp2[Nt][1][N];	//行列の積を計算するため(matrix_prod()でエラーが起こらないように)にr1[Nr][N]とr2[Nt][N]を3次元に変える
	double complex H_transpose[Nt][Nr][N];			//通信路行列Hのエルミート転置(Nt*Nr)
	double complex check[Nr][Nt][N];
	double complex W1[Nt][Nt][N];		//H^H×H (Nt×Nt)
	double complex W2[Nt][Nt][N];		//W1の逆行列 (Nt×Nt)
	double complex W[Nt][Nr][N];			//ZF等化の重みW
	
	for(j=0;j<Nr;j++){
		for(i=0;i<N;i++)
			tmp1[j][0][i]=r1[j][i];		//最初tmp1[j][1][i]と書いて爆死
	}
	
	///////////////////////////// W=H^-1 versinon  ///////////////////////////////////
	
	// inv_matrix(Nt,N,H,H_inv);		// 通信路行列の逆行列を求める　
	
	// matrix_prod(Nt,Nr,1,N,H_inv,tmp1,tmp2);		//Hの逆行列H^-1を左からr1に掛ける→r2(Nt×1)行列
	
	
	///////////////////////////// W=inv(H^H×H)×H^H /////////////////////////////
	Hermitian_transpose(Nr,Nt,N,H,H_transpose);		//　通信路行列Hのエルミート転置を求める 
	
	matrix_prod(Nt,Nr,Nt,N,H_transpose,H,W1);		//H^H(Nt×Nr)とH(Nr×Nt)の積W1(Nt×Nt)
	inv_matrix(Nt,N,W1,W2);							//W1の逆行列W2(Nt×Nt)
	matrix_prod(Nt,Nt,Nr,N,W2,H_transpose,W);		//W2(Nt×Nt)とH^H(Nt×Nr)の積W (Nt×Nr)
	matrix_prod(Nt,Nr,1,N,W,tmp1,tmp2);				//W(Nt×Nr)とtmp1(Nr×1)の積tmp2(Nt×1)
	
	
	///tmp2は3次元の配列なので2次元の配列に戻す//
	for(j=0;j<Nt;j++){
		for(i=0;i<N;i++)
			r2[j][i]=tmp2[j][0][i];		//tmp2[j][1][i]にしてはならない
	} 
	
	
	
	//確認　この確認でtmp[j][1][i]が間違いだと判明
	// for(i=0;i<N;i++){
		// for(j=0;j<Nr;j++){
			// printf("r1=%f\ttmp1=%f\n",creal(r1[j][i]),creal(tmp1[j][0][i]));
		// }
	// }
	
	
	//逆行列になっている確認
	// matrix_prod(Nt,Nr,Nt,N,H_inv,H,check);
	
	// putchar('\n');
	// for(i=0;i<N;i++){
		// for(k=0;k<Nt;k++){
			// for(j=0;j<Nr;j++){
				// printf("H=%f,H^-1=%f\n",creal(check[k][j][i]),creal(H_inv[k][j][i]));
			// }
		// }
	// }
	
	
	
}


////////////////////////    MMSE等化     ///////////////////////////////////////////
// 			　　　　　　		  通信路行列　　			受信側雑音電力		       等化前信号　　　				 等化後信号　　　
void MMSE(double complex H[Nr][Nt][N],double variance,double complex r1[Nr][N],double complex r2[Nt][N]){
	int i,j,k;
	double complex tmp1[Nr][1][N],tmp2[Nt][1][N];	//行列の積を計算するため(matrix_prod()でエラーが起こらないように)にr1[Nr][N]とr2[Nr][N]を3次元に変える
	double complex H_transpose[Nt][Nr][N];			//通信路行列H(Nr×Nt)のエルミート転置H^H(Nt×Nr)
	double complex W[Nt][Nr][N];		//MMSEフィルタの伝達関数

	for(j=0;j<Nr;j++){
		for(i=0;i<N;i++){
			tmp1[j][0][i]=r1[j][i];		//最初tmp1[j][1][i]と書いて爆死	
		}
	}
	
	//通信路行列Hのエルミート転置を求める
	Hermitian_transpose(Nr,Nt,N,H,H_transpose);		//H(Nr×Nt)のエルミート転置H(Nt×Nr)

	//方法①指導書 W=H^H*inv(H*H^H+σ^2*I_Nr)////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	/* double complex W1[Nr][Nr][N],W2[Nr][Nr][N];
	
	matrix_prod(Nr,Nt,Nr,N,H,H_transpose,W1);		//W1は(Nr×Nr)
	
	for(j=0;j<Nr;j++){
		for(i=0;i<Nr;i++){
			for(k=0;k<N;k++){//σを足すのは全成分ではなく，単位行列をかけてるから，対角行列の場所だけ
			// W1[i][j][k]+=variance;			//間違い！→対角成分の満たす必要がある,varianceは受信アンテナ全体の雑音電力だからNrいらない	
		
	// W1[i][j][k]+=variance*Nr;		//指導書通りならこれ間違い，sigmaが受信アンテナ全体の雑音電力を表すから
			if(i==j)
				W1[i][j][k]+=Nr*variance;		//i=jつまり対角成分にのみ雑音電力σ^2をたす，それ以外は足さない=成分は変化しない
			}
		}
	}
	
	inv_matrix(Nr,N,W1,W2);		//W2は(Nr×Nr)
	
	//MMSE線形フィルタの伝達関数Wを求める
	matrix_prod(Nt,Nr,Nr,N,H_transpose,W2,W);		//Wは(Nt×Nr)
	 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//方法②様々な文献にあるMMSE等化  W=inv(H^H*H+σ^2*I_Nr)*H^H  →inv()の中身の行列の形が異なるからミス////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* matrix_prod(Nt,Nr,Nt,N,H_transpose,H,W1);		//W1は(Nt×Nt)

	for(j=0;j<Nr;j++){
		for(i=0;i<Nr;i++){
			for(k=0;k<N;k++){//σを足すのは全成分ではなく，単位行列をかけてるから，対角行列の場所だけ
				if(i==j)
					W1[i][j][k]+=variance;		//i=jつまり対角成分にのみ雑音電力σ^2をたす，それ以外は足さない=成分は変化しない
			}
		}
	}
	
	inv_matrix(Nt,N,W1,W2);		//W2は(Nt×Nt)
	
	//MMSE線形フィルタの伝達関数Wを求める
	matrix_prod(Nt,Nt,Nr,N,W2,H_transpose,W);		//Wは(Nt×Nr)
	 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//方法③指導書＋参考文献 W=inv(H^H*H+σ^2*I_Nt)*H^H////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double complex W1[Nt][Nt][N],W2[Nt][Nt][N];
	
	matrix_prod(Nt,Nr,Nt,N,H_transpose,H,W1);		//W1は(Nt×Nt)
	
	for(j=0;j<Nt;j++){
		for(i=0;i<Nt;i++){
			for(k=0;k<N;k++){//σを足すのは全成分ではなく，単位行列をかけてるから，対角行列の場所だけ
			// W1[i][j][k]+=variance;			//間違い！→対角成分の満たす必要がある,varianceは受信アンテナ全体の雑音電力だからNrいらない	
			// W1[i][j][k]+=variance*Nr;		//指導書通りならこれ間違い，sigmaが受信アンテナ全体の雑音電力を表すから
			if(i==j)
				W1[i][j][k]+=variance;		//i=jつまり対角成分にのみ雑音電力σ^2をたす，それ以外は足さない=成分は変化しない
			}
		}
	}
	
	inv_matrix(Nt,N,W1,W2);		//W2は(Nt×Nt)
	
	//MMSE線形フィルタの伝達関数Wを求める
	matrix_prod(Nt,Nt,Nr,N,W2,H_transpose,W);		//Wは(Nt×Nr)
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//MMSE等化　部分
	matrix_prod(Nt,Nr,1,N,W,tmp1,tmp2);			//W(Nt×Nr)×tmp1(Nr×1)の積tmp2(Nt×1)
	
	//tmp2は3次元の配列なので2次元の配列に戻す
	for(j=0;j<Nt;j++){
		for(i=0;i<N;i++)
			r2[j][i]=tmp2[j][0][i];		//tmp2[j][1][i]にしてはならない
	} 

}

////////////////////////////////////////////順列を格納////////////////////////////////////////////////
int tmp[Q];
int k=0;

void copy_perm(int perm1[],int num[]){
	int i;

	for(i=0;i<Nt;i++){
		perm1[num[0]]=tmp[i];
		// printf("%d,perm=%d\n",k,perm1[k]);
		num[0]++;
	}

}

void perm_sub(int perm[],int n1,int n2,int n3[]){
	int i;
	
	if(n1==n2){
		copy_perm(perm,n3);
	}else{
		for(i=0;i<Q;i++){
			tmp[n2]=i;
			perm_sub(perm,n1,n2+1,n3);
		}
	}
	
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//           				(Nr×Nt)行列				 (Nr×1)行列
void MLD_QAM16(double complex H[Nr][Nt][N],double complex rt[Nr][N],double complex estimate[Nt][N]){
	int i,j,k,n,m;
	int n1[1];
	int p=(int)pow(Q,Nt)*Nt;
	int perm_num[p];				//数字の順列を格納 perm_num={0,0,0,1,0,2,0,3,1,0,･･･}
	int q=(int)pow(Q,Nt);			//q=Q^Nt
	double amp=1/sqrt(10*s);		//16QAMのIQ平面上の点の振幅
	double complex st_QAM16[Q]={amp+amp*I,-amp+amp*I,-amp-amp*I,amp-amp*I,3*amp+amp*I,-3*amp+amp*I,-3*amp-amp*I,3*amp-amp*I,amp+3*amp*I,-amp+3*amp*I,-amp-3*amp*I,amp-3*amp*I,3*amp+3*amp*I,-3*amp+3*amp*I,-3*amp-3*amp*I,3*amp-3*amp*I};	
	double complex replica[Nt][q][N];			//送信信号の全組み合わせ	縦にNt個，横にQ^Nt個ある行列
	double complex Hx[Nr][q][N];				//(Nr×Q^Nt)行列
	double min[N]={0};			
	double tmp;
	int key[N];
	
	
	
	//順列格納 (00)(01)(02)(03)...(0 15)(10)(11)(12)(13)...(1 15)(20)(21)(22)...(15 14)(15 15)
	n1[0]=0;
	perm_sub(perm_num,Nt,0,n1);		//perm_num={0,0,0,1,0,2,0,3,0,4,...0,15,1,0,1,1,1,2....15,15}
	
	
	// 確認
	// for(i=0;i<p;i++)
		// printf("%d\n",perm_num[i]);
	
	
	//jの値 012...Q^Nt
	//    ↓□□□
	//    ↓□□□　　jを固定して（縦に刻んで考えていく），すべての信号パターンを格納する

	for(k=0;k<N;k++){
		m=0;
		for(j=0;j<q;j++){
			for(i=0;i<Nt;i++){
				replica[i][j][k]=st_QAM16[perm_num[m]];
				m++;
			}
		}
	}
	
	// 通信路行列H(Nr×Nt)とreplica信号(Nt×Q^Nt)の積→出力の行列Hx(Nr×Q^Nt)
	matrix_prod(Nr,Nt,q,N,H,replica,Hx);
	
	//全探査
	for(j=0;j<N;j++){
		for(i=0;i<q;i++){		//入力信号の候補を順番に試す
			tmp=0.0;			//候補ごとにtmpを初期化する
			for(k=0;k<Nr;k++){
				if(i==0){
					min[j]+=cabs(rt[k][j]-Hx[k][i][j]);			//最初は初期値を入れる,y-Hxの結果Nr×1の行列ができるがNr個分の受信信号r(t)との距離を足す
					key[j]=0;		//最小値の時の候補の番号
				}
				else{			//i=1以降
					tmp+=cabs(rt[k][j]-Hx[k][i][j]);		//k=0～Nr個分距離を足す
				}
					
			}//受信アンテナ数Nrだけ距離を計算
			
			if(i>0){	//i=0の時はtmp=0.0でminが0.0に更新されてしまう
				if(min[j]>tmp){
						min[j]=tmp;		//最小値を更新する y1-Hx，y2-Hx，y3-Hx...と距離を足して最小値
						key[j]=i;		//最小値の時の候補番号を格納
				}
			}
		}//送信信号候補で回していく
	}//サブキャリア数Nだけ回す
	
	//推定値をコピー
	for(j=0;j<N;j++){
		for(i=0;i<Nt;i++)
			estimate[i][j]=replica[i][key[j]][j];		//候補信号の中から距離が最小となる番号keyをコピー
	}
	
	

	// printf("確認\n");
	// for(j=0;j<q;j++){
		// for(i=0;i<Nt;i++){
			// printf("%f\t%f\n",creal(Hx[i][j][0]),cimag(Hx[i][j][0]));
			// printf("%f\t%f\n",creal(replica[i][j][0]),cimag(replica[i][j][0]));
		// }
	// }

}


