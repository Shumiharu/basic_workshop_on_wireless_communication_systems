/*/////////////////////  LDPC符号を施した時のBER特性　　///////////////////////////////
・４ループ最小の時の検査行列から生成される４ループ最小生成行列Gを書き込む

*/

#include <stdio.h>
#include <stdlib.h>		//乱数使用時のプログラム
#include <time.h>		//time関数(乱数の種を設定するため)使用時のプログラム
#include <math.h>		//数学的関数使用時のプログラム
#include <complex.h>	//複素数使用時のプログラム



#define s 1						//オーバーサンプル比
// #define K 799					//送信QPSK複素信号のシンボル数，符号語(テールビットを含めて3200になるように設定)


#define M 128					//パリティビット長
#define N 256				//符号長
#define bit_length N-M			//情報ビット数
#define column_weight 3			//列重みj
#define row_weight 6		//行重みk


//LDPC符号用
#define lmax 20		//最大繰り返し数
#define choise 2	//1=検査行列Hを自分で作成しBERを出す，2=検査行列Hを先輩の見つけた性能の良いHをそのまま使いBERを出す

//BER特性用↓↓↓↓
#define bit 2 				//1シンボル当たりのビット数 BPSKの場合は1
#define BnTs 1   				
#define first 1000	//一旦雑にBER出す用のパラメータ
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
	
	//行列Xの確認
	// putchar('\n');
	// for(i=0;i<N;i++){
		// putchar('\n');
		// for(j=0;j<N+M;j++)
			// printf("%d ",X[i][j]);
		
	// }
	
	// array_XOR(N+M,X[1],X[4]);
	// row_swap(1,10,N,M+N,X);
	// putchar('\n');
	// printf("行入れ替え後のX\n");
	// for(i=0;i<N;i++){
		// putchar('\n');
		// for(j=0;j<N+M;j++)
			// printf("%d ",X[i][j]);	
	// }
	

/////////行列Xの変形→生成行列Gの作成///////////////////////////////////////////////////////////////////////////////////////////////////////////
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
		
		
		// putchar('\n');
		// printf("%d回目変換後の行列X\n",j);
		// for(k=0;k<N;k++){
			// putchar('\n');
			// for(l=0;l<N+M;l++)
				// printf("%d ",X[k][l]);
		// }
	
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
					// putchar('\n');
					// printf("%d回目row_swap変換後の行列X\n",j);
					// for(k=0;k<N;k++){
						// putchar('\n');
						// for(l=0;l<N+M;l++)
							// printf("%d ",X[k][l]);
					// }
					break;		//入れ替えたらこのforループを抜ける
				}//行入れ替え end
				
				if(i==N-1){	//1を持つ行が見つからなかったら→列入れ替え
					for(m=M+N-1; m>M-1; m--){		//最後の列からさかのぼって探査，列シフト
						if(X[j-M][m]==1){		//最後の列からさかのぼって探査し，行j-Mに要素1を持つ列を見つける
							column_swap(j,m,N,M+N,X);	//1を持つ列(m列目)探査後，列jと列mを入れ替える
							column_swap(j-M,m-M,M,N,H);	//検査行列Hの対応する列を入れ替える
							// putchar('\n');
							// printf("%d回目column_swap列変換後の行列X,%d列目と%d列目入れ替え\n",j,j,m);
							// for(k=0;k<N;k++){
								// putchar('\n');
								// for(l=0;l<N+M;l++)
									// printf("%d ",X[k][l]);
							// }
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
	
		// putchar('\n');
		// printf("%d回目変換後の行列X\n",j);
		// for(k=0;k<N;k++){
			// putchar('\n');
			// for(l=0;l<N+M;l++)
				// printf("%d ",X[k][l]);
		// }
	}//列シフト end (行列Xの右側変形終了)
	
//////(3)行列Xの一部から生成行列Gを抜き出す///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(i=M; i<N; i++){
		for(j=M; j<M+N; j++){
			G[i-M][j-M]=X[i][j];
		}
	}
	
	sprintf(fname, "4loopmin_Gmatrix,N=%d,M=%d,%d,%d.txt",N,M,column_weight,row_weight);
	fp=fopen(fname,"w");
	if(fp==NULL){
		printf("%sERROR!\n",fname);
		exit(EXIT_FAILURE);
	}
	for(i=0; i<N-M; i++){//生成行列Gの書き込み
		for(j=0; j<N; j++){
			fprintf(fp,"%d ", G[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	printf("%s Finish Writting\n",fname);

//////(4)GH^t=0の確認///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int check[N-M][M];	//check=GH^t=0が成り立っているか確認
	
	for(i=0; i<N-M; i++){
		for(j=0; j<M; j++)
			check[i][j]=0;
	}

	for(i=0; i<(N-M); i++){
		for(j=0; j<M; j++){
			for(k=0; k<N; k++){
				check[i][j]^=G[i][k]*H[j][k];
				
			}
		}
	}
	
	//GH^t=0の確認
	// putchar('\n');
	// printf("GH^t=0の確認\n");
	// for(i=0;i<N-M;i++){
		// putchar('\n');
		// for(j=0;j<M;j++)
			// printf("%d ",check[i][j]);	
	// }
	
	
	//変換後行列Xの確認
	// putchar('\n');
	// printf("変換後の行列X\n");
	// for(i=0;i<N;i++){
		// putchar('\n');
		// for(j=0;j<N+M;j++)
			// printf("%d ",X[i][j]);
	// }
	
	//列入れ替え後の検査行列Hの確認
	// putchar('\n');
	// printf("変換後の行列H\n");
	// for(i=0;i<M;i++){
		// putchar('\n');
		// for(j=0;j<N;j++)
			// printf("%d ",H[i][j]);	
	// }
	
	//列入れ替え後の検査行列Hの転置の確認
	// putchar('\n');
	// printf("変換後の行列Hの転置\n");
	// for(j=0;j<N;j++){
		// putchar('\n');
		// for(i=0;i<M;i++)
			// printf("%d ",H[i][j]);	
	// }
	
	//生成行列Gの確認
	// putchar('\n');
	// printf("生成行列G\n");
	// for(i=0;i<N-M;i++){
		// putchar('\n');
		// for(j=0;j<N;j++)
			// printf("%d ",G[i][j]);
		
	// }
	
	// int m_bit[N-M];		//情報ビットm（1×(N-M)）
	int (*m_bit)=(int(*))malloc((N-M)*sizeof(int));
	// int codeword[N];	//符号語c（1×N），符号長N
	int (*codeword)=(int(*))malloc(N*sizeof(int));
	// int check1[M];		//cH^t=0の確認用
	int (*check1)=(int(*))malloc(M*sizeof(int));
	// double complex BPSK_signal[N];		//ベースバンドBPSK信号
	double complex (*BPSK_signal)=(double complex(*))malloc(N*sizeof(double complex));
	// double complex QPSK_signal[N/2];	//ベースバンドQPSK信号
	double complex (*QPSK_signal)=(double complex(*))malloc((N/2)*sizeof(double complex));
	// double complex awgn[N/2];					//AWGN
	double complex (*awgn)=(double complex(*))malloc((N/2)*sizeof(double complex));
	// double complex receiving_signal[N/2];		//受信ベースバンドQPSK信号
	double complex (*receiving_signal)=(double complex(*))malloc((N/2)*sizeof(double complex));
	
	double sigma;				//ガウス雑音の分散			
	double EbN0;					//Eb/N0
	double count;				//BERカウント用
	double ber;  				//誤り率
	int loop,LOOP;
	
	
	// int N_m[M][row_weight];			//Nm，行mと固定したときのHmn=1を満たす行nの集合
	int (*N_m)[row_weight]=(int(*)[row_weight])malloc(M*row_weight*sizeof(int));
	// int M_n[N][column_weight];		//Mn，列nと固定したときのHmn=1を満たす列mの集合
	int (*M_n)[column_weight]=(int(*)[column_weight])malloc(N*column_weight*sizeof(int));
	// double LLR[N];		//受信LLR値
	double (*LLR)=(double(*))malloc(N*sizeof(double));
	// double zizenLLR[M][N][lmax];		//事前LLR umn
	double (*zizenLLR)[N][lmax]=(double(*)[N][lmax])malloc(M*N*lmax*sizeof(double));
	// double gaibuLLR[M][N][lmax];		//外部LLR vmn
	double (*gaibuLLR)[N][lmax]=(double(*)[N][lmax])malloc(M*N*lmax*sizeof(double));
	
	double tmp1;	//Π{sign(λ_n'+u_mn'(l)}用
	double tmp2;	//Σ{f(|λ_n'+u_mn'(l)|)}用
	double tmp3;	//Σ{v_m'n(l)}用
			
	// int parity_check[M];	//C_hat*H^t=
	int (*parity_check)=(int(*))malloc(M*sizeof(int));
	int sum;
	
	
	// int c_hat[N];		//一時推定語
	int (*c_hat)=(int(*))malloc(N*sizeof(int));
	// int estimate[N-M];		//推定情報ビット
	int (*estimate)=(int(*))malloc((N-M)*sizeof(int));
	double tmp4;		//Σ{v_mn(l)}用
						
	
	
	
//////↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓BER算出↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	printf("  「LDPC符号のBER特性」　　情報ビット:%d　，符号長：%d , 符号化率=%.2f  lmax:%d \n",N-M,N,(double)(N-M)/N,lmax);
	for(EbN0=minEbN0; EbN0<=maxEbN0; EbN0+=0.5){
		count=0.0;			//BERカウント初期化
		LOOP=first;			//Eb/N0の値ごとにfirst回でBERを軽く算出する
		//☆☆☆☆☆☆☆☆☆雑音の大きさ
		sigma=(double)BnTs/(2.0*s*bit)*pow(10.0,(-0.1*EbN0))*((double)N/(double)(N-M));		//ガウス雑音の分散算出，
		for(loop=0; loop<LOOP; loop++){	
			
			//////符号語cの生成///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			//////(1)情報ビットmの生成///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			
			for(i=0; i<N-M; i++)
				m_bit[i]=rand()%2;
			
			//情報ビットmの確認
			// putchar('\n');
			// putchar('\n');
			// printf("情報ビットmの確認\n");
			// for(i=0;i<N-M;i++)
				// printf("%d ",m_bit[i]);	
			
			//////(2)生成行列Gをかけて符号語c(符号長N)を作る，///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			//初期化
			for(i=0; i<N; i++)
				codeword[i]=0;
			

			for(j=0; j<N; j++){		//列シフト
				for(i=0; i<N-M; i++){	//行シフト
					codeword[j]^=m_bit[i]*G[i][j];		//c=mGにより符号語生成
				}
			}

			//符号語cの確認
			// putchar('\n');
			// putchar('\n');
			// printf("符号語cの確認(組織符号だから情報ビットがでるはず)\n");
			// for(i=0;i<N;i++)
				// printf("%d ",codeword[i]);	
			
			//cH^t=0の確認
		/* 	// for(i=0; i<M; i++)
				// check1[i]=0;
			
			// for(i=0; i<M; i++){
				// for(j=0; j<N; j++){
					// check1[i]^=codeword[j]*H[i][j];
				// }
			// }
					
			//cH^t=0の確認
			// putchar('\n');
			// printf("cH^t=0の確認\n");
			// for(i=0;i<M;i++)
				// printf("%d ",check1[i]);	
			 */
			
			//////BPSK QPSK modulation///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			//BPSK変調
			// for(i=0; i<N; i++){
				// if(codeword[i]==0)	BPSK_signal[i]=1.0;
				// if(codeword[i]==1)	BPSK_signal[i]=-1.0;
			// }
			
			//BPSKベースバンド信号の確認
			// putchar('\n');
			// putchar('\n');
			// printf("BPSKベースバンド信号の確認\n");
			// for(i=0;i<N;i++)
				// printf("%.1f ",creal(BPSK_signal[i]));	
			
			//QPSK変調
			bit_QPSK(N,codeword,QPSK_signal);

			//QPSKベースバンド信号の確認
			// putchar('\n');
			// putchar('\n');
			// printf("QPSKベースバンド信号の確認\n");
			// for(i=0;i<N/2;i++)
				// printf("%.1f %.1f\n",creal(QPSK_signal[i]),cimag(QPSK_signal[i]));	
			
			
			//////QPSKベースバンド信号がAWGN通信路通過///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
			//////（１）AWGN生成///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			boxmuller(N/2,awgn,0.0,sigma);
		
			//////（2）AWGN付加///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			for(i=0; i<N/2; i++){
				receiving_signal[i]=QPSK_signal[i]+awgn[i];
				// receiving_signal[i]=QPSK_signal[i];		//AWGNなし
			}
			
			//受信されたQPSKベースバンド信号の確認
			// putchar('\n');
			// putchar('\n');
			// printf("受信されたQPSKベースバンド信号の確認\n");
			// for(i=0;i<N/2;i++)
				// printf("%f %f \n",creal(receiving_signal[i]),cimag(receiving_signal[i]));	
			
			//////sum-product復号の準備///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
			// int N_m[M][row_weight];			//Nm，行mと固定したときのHmn=1を満たす行nの集合
			// int M_n[N][column_weight];		//Mn，列nと固定したときのHmn=1を満たす列mの集合
			// double LLR[N];		//受信LLR値
		
			//////（1）N_m生成///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			for(i=0; i<M; i++){	//行シフト
				k=0;
				for(j=0; j<N; j++){	//列シフト
					if(H[i][j]==1){	
						N_m[i][k]=j;
						k++;
					}
				}
			}
				
			//N_mの確認
			// putchar('\n');
			// putchar('\n');
			// printf("N_m(行mを固定したときのHmnを満たす列nの集合)の確認\n");
			// for(i=0;i<M;i++){
				// putchar('\n');
				// for(j=0; j<row_weight; j++)
					// printf("%d ",N_m[i][j]);	
			// }
			
			//////（2）M_n生成///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			for(i=0; i<N; i++){	//列シフト
				k=0;
				for(j=0; j<M; j++){	//行シフト
					if(H[j][i]==1){	
						M_n[i][k]=j;
						k++;
					}
				}
			}			
					
			//M_nの確認
			// putchar('\n');
			// putchar('\n');
			// printf("M_n(列nを固定したときのHmnを満たす行mの集合)の確認\n");
			// for(i=0;i<N;i++){
				// putchar('\n');
				// for(j=0; j<column_weight; j++)
					// printf("%d ",M_n[i][j]);	
			// }
			
			//////（3)受信LLR算出///////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			for(i=0; i<N/2; i++)
				receiving_signal[i]/=sqrt(2.0);		//正規化
			
			for(i=0; i<N/2; i++){
				LLR[2*i]=2.0*creal(receiving_signal[i])/sigma;
				LLR[2*i+1]=2.0*cimag(receiving_signal[i])/sigma;
			}
			
			//受信LLR値確認
			// putchar('\n');
			// putchar('\n');
			// printf("受信LLR値確認\n");
			// for(i=0;i<N;i++){
				// printf("%d : %f\n",i,LLR[i]);	
			// }
			
			//////sum-product復号///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
			// double zizenLLR[M][N][lmax];		//事前LLR umn
			// double gaibuLLR[M][N][lmax];		//外部LLR vmn
			
			//////(1)事前LLRと外部LLRの初期化///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
			for(i=0; i<M; i++){
				for(j=0; j<N; j++){
					zizenLLR[i][j][0]=0.0;
					gaibuLLR[i][j][0]=0.0;
				}
			}	//l=0　end	
			
			// double tmp1;	//Π{sign(λ_n'+u_mn'(l)}用
			// double tmp2;	//Σ{f(|λ_n'+u_mn'(l)|)}用
			// double tmp3;	//Σ{v_m'n(l)}用
			
			////////////☆☆☆☆☆☆l<=lmax?
			for(l=1; l<=lmax; l++){
				//////(2)行処理(チェックノード処理）///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
				//i→m，j→n，N_m[i][k]→n'に対応する
				for(i=0; i<M; i++){	//チェックノードiにおいて
					for(j=0; j<row_weight; j++){	//N_m[i][j]：チェックノードiにおけるj番目に接続されている変数ノード
						tmp1=1.0;
						tmp2=0.0;
						for(k=0; k<row_weight;k++){
							if(N_m[i][k]!=N_m[i][j]){	//今注目してる変数ノードN_m[i][j]以外の変数ノードが演算に参加する
								tmp1*=sign(LLR[N_m[i][k]]+zizenLLR[i][N_m[i][k]][l-1]);
								tmp2+=Gallager_function(fabs(LLR[N_m[i][k]]+zizenLLR[i][N_m[i][k]][l-1]));
							}
						}//n'の分の計算 end
						// if(i==0&&j==0){	printf("tmp1=%f  tmp2=%f\n",tmp1,tmp2);}
						gaibuLLR[i][N_m[i][j]][l]=tmp1*Gallager_function(tmp2);
					}//変数ノードjの外部LLR値計算 end
				}//チェックノードiにおけるすべての外部LLR値計算 end
				
				//外部LLR値の確認
				// putchar('\n');
				// putchar('\n');
				// printf("外部LLR値の確認\n");
				// for(i=0; i<M; i++){
					// putchar('\n');
					// for(j=0; j<row_weight; j++)
					// printf("%d : %d : %f\n",i,N_m[i][j],gaibuLLR[i][N_m[i][j]][1]);
				// }
				
				//////(3)列処理(変数ノード処理）///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
				//i→n，j→m，M_n[i][k]→m'に対応する
				for(i=0; i<N; i++){	//変数ノードiにおいて
					for(j=0; j<column_weight; j++){	//M_n[i][j]：変数ノードiにおけるj番目に接続されているチェックノード
						tmp3=0.0;
						for(k=0; k<column_weight;k++){
							if(M_n[i][k]!=M_n[i][j]){	//今注目してるチェックノードM_n[i][j]以外のチェックノードが演算に参加する
								tmp3+=gaibuLLR[M_n[i][k]][i][l];
							}
						}//m'の分の計算 end
						// if(i==0&&j==0){	printf("tmp3=%f\n",tmp3);}
						zizenLLR[M_n[i][j]][i][l]=tmp3;
					}//変数ノードjの外部LLR値計算 end
				}//チェックノードiにおけるすべての外部LLR値計算 end
				
				//事前LLR値の確認
				// putchar('\n');
				// putchar('\n');
				// printf("事前LLR値の確認\n");
				// for(i=0; i<N; i++){
					// putchar('\n');
					// for(j=0; j<column_weight; j++)
						// printf("%d : %d : %f\n",M_n[i][j],i,zizenLLR[M_n[i][j]][i][1]);
				// }
				
				//////(4)一次推定語算出///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
				// int c_hat[N];		//一時推定語
				// int estimate[N-M];		//推定情報ビット
				// double tmp4;		//Σ{v_mn(l)}用
				
				for(i=0; i<N; i++){	//変数ノードi シフト
					tmp4=0.0;
					for(j=0; j<column_weight; j++){		//変数ノードiと接続しているチェックノードM_n[i][j]　シフト
						tmp4+=gaibuLLR[M_n[i][j]][i][l];
					}
					if(sign(LLR[i]+tmp4)>0)	c_hat[i]=0;
					if(sign(LLR[i]+tmp4)<0)	c_hat[i]=1;
				}
				
				//一時推定語c_hatの後半部分の情報ビットを抜き出す
				for(i=0; i<N-M; i++)
					estimate[i]=c_hat[i+M];
				
				
				//一時推定語c_hatの確認
				// putchar('\n');
				// putchar('\n');
				// printf("一時推定語c_hatの確認\n");
				// for(i=0;i<N;i++)
					// printf("%d ",c_hat[i]);	
				
				//推定情報ビットの確認
				// putchar('\n');
				// putchar('\n');
				// printf("推定情報ビットの確認\n");
				// for(i=0;i<N-M;i++)
					// printf("%d ",estimate[i]);	
				
				
				
				//////(5)パリティ検査///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
				// int parity_check[M];	//C_hat*H^t=0の確認
				// int sum;
				
				//初期化
				for(i=0; i<M; i++)
					parity_check[i]=0;
				
				//パリティ検査c_hat*H^t
				for(i=0; i<M; i++){
					for(j=0; j<N; j++){
						parity_check[i]^=c_hat[j]*H[i][j];
					}
				}
				
				//アルゴリズム終了かどうかの判定(parity_checkの全成分を足して(XORじゃない普通の足し算)，総和が0なら終了)
				sum=0;
				for(i=0; i<M; i++)
					sum+=parity_check[i];	//C_hat*H^tの計算により出てきた成分の和
				
				// printf("sum=%d\n",sum);
				if(sum==0)	break;	//全成分0ならsum-product endのループをbreak
				
				// cH^t=0の確認
				// putchar('\n');
				// printf("c_hat*H^t=0の確認\n");
				// for(i=0;i<M;i++)
					// printf("%d ",parity_check[i]);	
			 
				//////(6)最大繰り返し数lmaxに到達時///////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
				//for(l=1; l<lmax; l++)のループでl=lmaxになったときはアルゴリズムが終了するようになっている
			
				
				
				
			}//sum-product end
			
			
			//判定
			for(i=0; i<N-M; i++){
				if(m_bit[i]!=estimate[i])
					count++;
			}
			
			
			//BERの値によってループ調整
			if(loop==LOOP-1){
					ber=count/((double)loop*(N-M));		//一旦BERを算出
				
				if(ber>=1E-1 && ber<1)
					LOOP=pow(10,3);			//loop=first=5000→このifに入る→LOOP=10が代入される→最初のforの部分でloop++される→loop<LOOPの判定でloopの方が大きいのでこのforを抜け出せる→BERを計算して結果表示
				if(ber>=1E-2 && ber<1E-1)
					LOOP=pow(10,3);
				if(ber>=1E-3 && ber<1E-2)
					LOOP=pow(10,4);
				if(ber>=1E-4 && ber<1E-3)
					LOOP=5*pow(10,4);
				if(ber>=1E-5 && ber<1E-4)
					LOOP=pow(10,5);
				if(ber>=1E-6 && ber<1E-5)
					LOOP=5*pow(10,5);
				if(ber>=1E-7 && ber<1E-6)
					LOOP=pow(10,6);
				if(ber>=1E-8 && ber<1E-7)
					LOOP=5*pow(10,7);
				if(ber>=0 && ber<1E-8)
					LOOP=pow(10,8);
				// printf("BER=%e, loop=%d, LOOP=%d, Eb/N0=%.2f dB\n",ber,loop,LOOP,EbN0);
			} 
		ber=count/((double)(1+loop)*(N-M));		//☆☆☆☆☆☆
		printf("loop=%d         Eb/N0=%.2f         %e\r",loop,EbN0,ber);	
			
		}//試行loop end
		
		//☆☆☆☆☆☆BERの中のN
		ber=count/((double)loop*(N-M));		//☆☆☆☆☆☆
		printf("loop=%d         Eb/N0=%.2f         %e\n",loop,EbN0,ber);	
	
	
	}//Eb/N0ループ end
	
	//メモリ開放
	free(H);		free(random_column);		free(G);	
	free(X);		free(m_bit);				free(codeword);
	free(check1);	free(BPSK_signal);			free(QPSK_signal);
	free(awgn);		free(receiving_signal);		free(N_m);
	free(M_n);		free(LLR);					free(zizenLLR);
	free(gaibuLLR);	free(parity_check);			free(c_hat);
	free(estimate);
	
	
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














