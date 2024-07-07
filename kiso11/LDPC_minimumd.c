/*/////////////////////  LDPC符号の符号間最小距離算出　　///////////////////////////////
・first×first個の符号語を比較し，最小距離d_minを算出
・ｃ１，ｃ２が同じ符号語だと最小距離が０になってしまうからその時はdminの更新は行わない
*/

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム



#define s 1						//オーバーサンプル比
// #define K 799					//送信QPSK複素信号のシンボル数，符号語(テールビットを含めて3200になるように設定)


#define M 150					//パリティビット長
#define N 300				//符号長
#define bit_length N-M			//情報ビット数
#define column_weight 2			//列重みj
#define row_weight 4		//行重みk


//LDPC符号用
#define lmax 40		//最大繰り返し数
#define choise 2	//1=検査行列Hを自分で作成しBERを出す，2=検査行列Hを先輩の見つけた性能の良いHをそのまま使いBERを出す

//BER特性用↓↓↓↓
#define bit 2 				//1シンボル当たりのビット数 BPSKの場合は1
#define BnTs 1   				
#define first 100000	//一旦雑にBER出す用のパラメータ
#define minEbN0 0.0
#define maxEbN0 5.0


////////////////////////関数プロトタイプ宣言/////////////////////////////////
void QPSK(int nd, double complex base[]);
void inter(double complex s1[],double complex s2[]);		//先頭サンプルのみ抜き出す
void boxmuller(int num,double complex noise[],double mean,double variance);
void QPSK_bit(int nd,double complex input[],int output[]);
void bit_QPSK(int nd, int input[], double complex output[]);
int arg_min(int n1, int n2);
void matrix_prod(int n1,int n2,int n3,double complex mx1[n1][n2],double complex mx2[n2][n3],double complex mxout[n1][n3]);
void inv_matrix(int num,int nd,double complex IN[num][num][nd],double complex OUT[num][num][nd]);

/////////// 順列(permutation)0～num-1までのnum個の数字をランダムに入れ替える関数 /////////////
void shuffle(int num, int perm[]){
	int i,j;
	int tmp;
	
	//順列に0～num-1の整数を格納
	for(i=0; i<num; i++)
		perm[i]=i;
	
	//順列の要素，整数を並び替える
	for(i=0; i<num; i++){
		j=rand()%num;		//余り，0～numまでの値がjに格納される
		tmp=perm[i];		//一時的にコピー
		perm[i]=perm[j];	//要素iと要素jを入れ替える
		perm[j]=tmp;		//要素jにtmpを入れる
	}
}

////////////////　nr×nc行列matrixにおいて，指定した2つの行n1と行n2を入れ替える関数 ////////////////
void row_swap(int n1,int n2,int nr,int nc,int matrix[nr][nc]){
	int i,j,k;
	int tmp[nc];		//行列matrixの指定行n1の要素をコピー
	
	//行n1をtmpにコピー
	for(i=0; i<nc; i++){
		tmp[i]=matrix[n1][i];
	}
	
	for(i=0; i<nc; i++){
		matrix[n1][i]=matrix[n2][i];
		matrix[n2][i]=tmp[i];
	}
}

////////////////　nr×nc行列matrixにおいて，指定した2つの列n1と列n2を入れ替える関数 ////////////////
void column_swap(int n1,int n2,int nr,int nc,int matrix[nr][nc]){
	int i,j,k;
	int tmp[nr];		//行列matrixの指定列n1の要素をコピー
	
	//列n1をtmpにコピー
	for(i=0; i<nr; i++){
		tmp[i]=matrix[i][n1];
	}
	
	for(i=0; i<nr; i++){
		matrix[i][n1]=matrix[i][n2];
		matrix[i][n2]=tmp[i];
	}
}

////////////////　行数または列数numの行列において，行または列2つの排他的論理和をとる関数 ////////////////
//			num:行数または列数    data1をdata2に加える(=XOR)
void array_XOR(int num, int data1[num], int data2[num]){
	int i,j;
	
	for(i=0; i<num; i++)
		data2[i]=data2[i]^data1[i];
	
}

double Gallager_function(double x){
	double fx;
	
	if(x<0.0000001){
		fx=(log((exp(0.0000001)+1.0)/(exp(0.0000001)-1.0)));
	}else if(x>30.0){
		fx=(log((exp(30.0)+1.0)/(exp(30.0)-1.0)));
	}else{
		fx=(log((exp(x)+1.0)/(exp(x)-1.0)));
	}
	return fx;
}

double sign(double x){
	
	if(x<0)
		return -1.0;
	else
		return 1.0;	
	
}



int main(void)
{
	
	
	
///////////////////////      乱数の種を設定     ///////////////
	srand((unsigned int)time(NULL)); 
	
	int i,j,k,l,m,n,o;
	

////////choise=1 : 検査行列Hの生成 //////////////////////////////////////////////////////////////////////////////////////////////////////
	// int H[M][N];		//検査行列H(M×N)
	int (*H)[N]=(int(*)[N])malloc(M*N*sizeof(int));		//検査行列H(M×N)
	// int random_column[N];		//検査行列Hの2段目以降の要素をランダムに入れ替えるため
	int (*random_column)=(int(*))malloc(N*sizeof(int));
	
	
	if(choise==1){
		//検査行列Hの1段目のサブブロックに階段状に1を格納
		k=0;
		for(l=0;l<N/row_weight;l++){	//N/row_weight回
			for(i=row_weight*l;i<row_weight*(l+1);i++){	//行重みkごとに同じ操作を繰り替えす
				for(j=0;j<M/column_weight;j++){		//1段目のサブブロック作成，縦向きに順に作成していく
					H[j][i]=(j==k)? 1:0;		//j=k(k段目にのみ1)，それ以外は0を格納
				}
			}
			k++;		//1を階段状に入れるためのk++
		}
		
		//検査行列Hの2段目以降のサブブロックの格納
		for(i=1; i<column_weight; i++){
			shuffle(N,random_column);		//列をランダムに入れ替えるため
			for(j=0; j<N; j++){
				for(k=0; k<M/column_weight; k++){
					H[(M/column_weight)*i+k][j]=H[k][random_column[j]];
				}//サブブロックの行数kループ　end
			}//検査行列Hの列数j end
		}//1段目以降のサブブロックiループ end
		
		//行入れ替えのためのランダム列の確認
		// putchar('\n');
		// for(i=0; i<N; i++)
			// printf("%d ",random_column[i]);	
		
		//検査行列Hの確認
		// putchar('\n');
		// for(i=0;i<M;i++){
			// putchar('\n');
			// for(j=0;j<N;j++)
				// printf("%d ",H[i][j]);
		// }
		// putchar('\n');
					
		
		// FILE *fp;//ファイル作成
		// char fname[100];
		// sprintf(fname, "Hmatrix,N=%d,M=%d,%d,%d.txt",N,M,column_weight,row_weight);
		// fp=fopen(fname,"w");
		// if(fp==NULL){
			// printf("%sERROR!\n",fname);
			// exit(EXIT_FAILURE);
		// }
		// for(i=0; i<M; i++){//検査行列Hの書き込み
			// for(j=0; j<N; j++){
				// fprintf(fp,"%d ", H[i][j]);
			// }
			// fprintf(fp,"\n");
		// }
		// fclose(fp);
		// printf("%s Finish Writting\n",fname);

	}//choise=1 end
	
///////////choise=2 : 検査行列Hの読み込み(先輩の見つけた性能の良い検査行列を読み込み，その検査行列でBERを出す)/////////////////////////////////////////////////////////////////////////////////////////////////////	
	FILE *fp;//ファイル読み込み
	char fname[100];
	
	if(choise==2){
		// sprintf(fname, "Hmatrix2,N=%d,M=%d,j=%d,k=%d.txt",N,M,column_weight,row_weight);
		// sprintf(fname, "bestcheckmatrix,N=%d,M=%d,colweight=%d,rowweight=%d.txt",N,M,column_weight,row_weight);		//奥村先輩の最適な検査行列
		sprintf(fname, "4loopmin_Hmatrix,N=%d,M=%d,j=%d,k=%d.txt",N,M,column_weight,row_weight);		//自分の4ループ最小の検査行列
		// char fname[] = "Hmatrix,N=600,M=300,3,6.txt";
		fp = fopen(fname,"r");		//読み込み(read)の"r"		//foepnでfnameと同じファイルを探す
		if(fp==NULL){
			printf("%s file not open!\n",fname);
			return -1;
		}else{
			printf("%s file opened!\n",fname);
		}
		
		
		for(i=0;i<M;i++){	
			for(j=0;j<N;j++){
				fscanf(fp,"%d",&H[i][j]);
			}
		}
		fclose(fp);
		// for(i=0; i<M; i++){//検査行列Hの表示
			// for(j=0; j<N; j++){
				// printf("%d ", H[i][j]);
			// }
			// putchar('\n');
		// }
		// putchar('\n');
	}//choise=2 end
	
	
////////////////行列Xの作成///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// int X[N][M+N];		//X=[H^t(N×M),I_N(N×N)]
	int (*X)[M+N]=(int(*)[M+N])malloc(N*(M+N)*sizeof(int));
	///Xの左側，検査行列Hの転置を格納
	for(i=0; i<N; i++){	//行数N
		for(j=0; j<M; j++)	//列数M
			X[i][j]=H[j][i];		//検査行列Hの転置を格納
	}
	
	//Xの右側，単位行列I_Nを格納
	for(i=0; i<N; i++){	//行数N
		for(j=M; j<M+N; j++){	//列数M～M+N-1
			if(i==j-M)	X[i][j]=1;
			else		X[i][j]=0;
		}
	}

	

/////////行列Xの変形→生成行列Gの作成///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// int G[N-M][N];		//生成行列G((N-M)×N)
	int (*G)[N]=(int(*)[N])malloc((N-M)*N*sizeof(int));
//////(1)行列Xの左側の変形///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(j=0; j<M; j++){	//列シフトj
		
		if(X[j][j]==1){		//対角成分が1ならば
			for(i=0; i<N; i++){	//行シフトi
				if(j!=i && X[i][j]==1){	//今見てる対角成分以外の行かつその成分X[j][i]=1のとき
					array_XOR(M+N,X[j],X[i]);	//対角成分をたすことで，その列jの1を消す
				}
			}//行シフト＆1消去 end
		}
		else if(X[j][j]==0){
			// printf("対角成分0\n");
			for(i=j+1; i<N; i++){	//今見てる行(=j)より下の行(j+1以降)
				if(X[i][j]==1){		//今見てる行より下の行の中で，1を持つ行を探査
					row_swap(j,i,N,M+N,X);	//1を持つ行を探査後，行を入れ替える
					// putchar('\n');
					// printf("%d回目swap変換後の行列X\n",j);
					// for(k=0;k<N;k++){
						// putchar('\n');
						// for(l=0;l<N+M;l++)
							// printf("%d ",X[k][l]);
			
					// }
					break;		//入れ替えたらこのforループを抜ける
				}//行入れ替え end
				if(i==N-1){	//1を持つ行が見つからなかったら→列入れ替え
					for(m=M+N-1; m>j; m--){		//最後の列からさかのぼって探査，列シフト
						if(X[j][m]==1){		//最後の列からさかのぼって探査し，行jに要素1を持つ列を見つける
							column_swap(j,m,N,M+N,X);	//1を持つ列(m列目)探査後，列jと列mを入れ替える
							break;
						}
					}//列を入れ替えたらbreakする
				}//列入れ替え end
			}//行探査＆入れ替え or 列探査＆入れ替え end
			
			///今見ている列jの中の1を消去するフェーズ
			for(i=0; i<N; i++){	//行シフトi
				if(j!=i && X[i][j]==1){	//今見てる対角成分以外の行かつその成分X[j][i]=1のとき
					array_XOR(M+N,X[j],X[i]);	//対角成分をたすことで，その列jの1を消す
				}
			}//行シフト＆1消去 end
		}
	
	}//列シフト end (行列Xの左側変形終了)
	
//////(2)行列Xの右側の変形///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(j=2*M; j<N+M; j++){	//列シフトj ☆この置き方だと行はj-M，列はjと表される☆
		if(X[j-M][j]==1){		//右下のスタートが，(行=M，列2M)，その対角成分X[j-M][j]が1ならば
			for(i=0; i<N; i++){	//行シフトi
				if((j-M)!=i && X[i][j]==1){	//今見てる対角成分以外の行かつその成分X[j][i]=1のとき
					array_XOR(M+N,X[j-M],X[i]);	//対角成分をたすことで，その列jの1を消す
				}
			}//行シフト＆1消去 end
		}
		else if(X[j-M][j]==0){		//今見ている対角成分X[j-M][j]=0なら
			// for(i=j-M+1; i<N; i++){	//今見てる行(=j-M)より下の行(j-M+1以降)
			for(i=j-M; i<N; i++){	//今見てる行(=j-M)以下の行(j-M+1以降)
				if(X[i][j]==1){		//今見てる行より下の行の中で，1を持つ行を探査できたなら
					row_swap(j-M,i,N,M+N,X);	//1を持つ行を探査後，行を入れ替える

					break;		//入れ替えたらこのforループを抜ける
				}//行入れ替え end
				
				if(i==N-1){	//1を持つ行が見つからなかったら→列入れ替え
					for(m=M+N-1; m>M-1; m--){		//最後の列からさかのぼって探査，列シフト
						if(X[j-M][m]==1){		//最後の列からさかのぼって探査し，行j-Mに要素1を持つ列を見つける
							column_swap(j,m,N,M+N,X);	//1を持つ列(m列目)探査後，列jと列mを入れ替える
							column_swap(j-M,m-M,M,N,H);	//検査行列Hの対応する列を入れ替える
							break;
						}
					}//列を入れ替えたらbreakする
				}//列入れ替え end
			}//行探査＆入れ替え or 列探査＆入れ替え end
			
			///今見ている列jの中の1を消去するフェーズ
			for(i=0; i<N; i++){	//行シフトi
				if((j-M)!=i && X[i][j]==1){	//今見てる対角成分以外の行かつその成分X[j][i]=1のとき
					array_XOR(M+N,X[j-M],X[i]);	//対角成分をたすことで，その列jの1を消す
				}
			}//行シフト＆1消去 end
		}
	

	}//列シフト end (行列Xの右側変形終了)
	
//////(3)行列Xの一部から生成行列Gを抜き出す///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(i=M; i<N; i++){
		for(j=M; j<M+N; j++){
			G[i-M][j-M]=X[i][j];
		}
	}
	

	
	int (*m_bit1)[first]=(int(*)[first])malloc((N-M)*first*sizeof(int));	//情報ビット1
	int (*m_bit2)[first]=(int(*)[first])malloc((N-M)*first*sizeof(int));	//情報ビット1
	int (*codeword1)[first]=(int(*)[first])malloc(N*first*sizeof(int));		//符号語1
	int (*codeword2)[first]=(int(*)[first])malloc(N*first*sizeof(int));		//符号語1
	
	
	double sigma;				//ガウス雑音の分散			
	double EbN0;					//Eb/N0
	double count=0;				//BERカウント用
	double ber;  				//誤り率
	int loop,LOOP;
	int d_min;
	int tmp_d=0;
	
	
	
//////↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓BER算出↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	printf("  「LDPC符号のBER特性」　　情報ビット:%d　，符号長：%d , 符号化率=%.2f  lmax:%d \n",N-M,N,(double)(N-M)/N,lmax);
		
		
	//////符号語cの生成///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//////(1)first個の情報ビットm1,m2の生成///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	for(j=0; j<first; j++){
		for(i=0; i<N-M; i++){
			m_bit1[i][j]=rand()%2;
			m_bit2[i][j]=rand()%2;
		}
	}
	
	//情報ビットmの確認
	// putchar('\n');
	// putchar('\n');
	// printf("情報ビットm1の確認\n");
	// for(j=0; j<first; j++){
		// for(i=0;i<N-M;i++){
			// printf("%d ",m_bit1[i][j]);	
		// }
	// }
	
	//////(2)生成行列Gをかけて符号語c(符号長N)を作る，///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//初期化
	for(j=0; j<first; j++){
		for(i=0; i<N; i++){
			codeword1[i][j]=0;
			codeword2[i][j]=0;
		}
	}
	
	//符号語c1,c2作成
	for(k=0; k<first; k++){
		for(j=0; j<N; j++){		//列シフト
			for(i=0; i<N-M; i++){	//行シフト
				codeword1[j][k]^=m_bit1[i][k]*G[i][j];		//c=mGにより符号語生成
				codeword2[j][k]^=m_bit2[i][k]*G[i][j];		//c=mGにより符号語生成
			}
		}
	}

	//符号語の最小距離探査
	//符号語c1,c2において，c1を固定しc2をfirst回シフトしてその都度最小距離を算出
	for(i=0; i<first; i++){
		for(j=0; j<first; j++){
			tmp_d=0;
			for(k=0; k<N; k++){
				tmp_d+=codeword1[k][i]^codeword2[k][j];	
			}//符号語c1,c2のハミング距離探査 end
			
			if(tmp_d!=0){	//tmp_d=0だとc1,c2が同じ符号語だからd_minの更新はなし
				//d_minの初期化
				if(count==0){
					d_min=tmp_d;
					count++;
				}
				
				//d_minの更新
				if(d_min>tmp_d)
					d_min=tmp_d;
			}
		}
	}

	printf("符号の最小距離\n");
	printf("  %d\n",d_min);
	
	
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
void QPSK(int nd, double complex base[]){
		int i,j;
		double tmp_I;	//同相成分Iを一時的に格納
		double tmp_Q;	//直交成分Qを一時的に格納
		
		for(i=0; i<nd; i++){
			tmp_I=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1	
			tmp_Q=1/(sqrt(2*s))*(-2*(rand()%2)+1);	//電圧の大きさが1/√2s , 乱数が0なら１で1ならー1
			
			for(j=0; j<s; j++)
				base[s*i+j]=tmp_I+tmp_Q*I;				//複素ベースバンド信号u(t)=I(t)+jQ(t)
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


/////////////////// QPSK複素信号をビット系列を変換する関数 //////////////////
//             入力のデータ数     QPSK複素信号(nd)   ビット系列(nd*bit)
void QPSK_bit(int nd,double complex input[],int output[]){
	int i,j;
	
	for(i=0;i<nd;i++){
		if(creal(input[i])>0 && cimag(input[i])>0){		//第1象限(00)
			output[i*2]=0;
			output[i*2+1]=0;
		}
		if(creal(input[i])<0 && cimag(input[i])>0){		//第2象限(10)
			output[i*2]=1;
			output[i*2+1]=0;
		}
		if(creal(input[i])<0 && cimag(input[i])<0){		//第3象限(11)
			output[i*2]=1;
			output[i*2+1]=1;
		}
		if(creal(input[i])>0 && cimag(input[i])<0){		//第4象限(01)
			output[i*2]=0;
			output[i*2+1]=1;
		}	
	}
}

/////////////////////　ビット列をQPSK信号に変換する関数  //////////////////////
//			    ビット長   ビット系列(nd)        QPSK信号(nd/2) 
void bit_QPSK(int nd, int input[], double complex output[]){
	int i,j;
	
	for(i=0;i<nd/2; i++){
		if(input[i*2]==0 && input[i*2+1]==0)
			output[i]=1/(sqrt(2*s))+I*1/(sqrt(2*s));	//第1象限(00)=A+jA
		
		if(input[i*2]==1 && input[i*2+1]==0)
			output[i]=-1/(sqrt(2*s))+I*1/(sqrt(2*s));	//第2象限(10)=-A+jA
		
		if(input[i*2]==1 && input[i*2+1]==1)
			output[i]=-1/(sqrt(2*s))-I*1/(sqrt(2*s));	//第1象限(11)=-A-jA
		
		if(input[i*2]==0 && input[i*2+1]==1)
			output[i]=1/(sqrt(2*s))-I*1/(sqrt(2*s));	//第1象限(01)=A-jA
	}
}


int arg_min(int n1, int n2){
	if(n1<=n2)
		return n1;		
	else
		return n2;
}
	

/////////////////  num1×num2とnum2×num3行列の掛け算 出力はnum1×num3の行列 ///////////////////////////////// 
void matrix_prod(int n1,int n2,int n3,double complex mx1[n1][n2],double complex mx2[n2][n3],double complex mxout[n1][n3]){
	int i,j,k,l;
	
	
	for(j=0;j<n1;j++){
		for(k=0;k<n3;k++){
			mxout[j][k]=0.0+0.0*I;
		}
	}

	
	for(j=0;j<n1;j++){
		for(k=0;k<n3;k++){
			for(l=0;l<n2;l++){
				mxout[j][k]+=mx1[j][l]*mx2[l][k];
			}
		}
	}

}

//////////////////////////////////// num×numの行列の逆行列を求める関数 ////////////////////////////
//                       nd:データ数
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














