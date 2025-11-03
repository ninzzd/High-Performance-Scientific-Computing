#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include <math.h>
#include <time.h>

#define pi 3.14159265359
void quad1 (double *p, double *w,int ngl)
{
    if(ngl==1)
    {
        p[ngl-1]=0;w[ngl-1]=2;
    }
    if(ngl==2)
    {
        p[0]=-0.577350269189625765;
        p[1]=-p[0];
        w[0]=1.0;
        w[1]=w[0];
    }
    if(ngl==3)
    {
        p[0]=0.774596669241483;
        p[1]=0;
        w[0]=0.5555555555555556;
        p[2]=-p[0];
        w[1]=0.8888888888888889;
        w[2]=w[1];
    }
    if (ngl==4)
    {
        p[0]=-0.861136311594053;
        p[1]=-0.339981043584856;
        p[2]=-p[1];
        p[3]=-p[0];
        w[0]=0.347854845137454;
        w[1]=0.652145154862546;
        w[2]=w[1];
        w[3]=w[0];

    }
}


void shapefn (double rvalue, double svalue, double *shape, double * dhdr, double *dhds)
{
    //  printf("r %lf   s %lf\n",rvalue,svalue );
    shape[0]=0.25*(1-rvalue)*(1-svalue);
    shape[1]=0.25*(1+rvalue)*(1-svalue);
    shape[2]=0.25*(1+rvalue)*(1+svalue);
    shape[3]=0.25*(1-rvalue)*(1+svalue);
    int i,j;

    //  printf("shape  %lf\n",shape[i] );

    dhdr[0]=  -0.25*(1-svalue);
    dhdr[1]=  0.25*(1-svalue);
    dhdr[2]=   0.25*(1+svalue);
    dhdr[3]=   -0.25*(1+svalue);

    dhds[0]=-0.25*(1-rvalue);
    dhds[1]=-0.25*(1+rvalue);
    dhds[2]=0.25*(1+rvalue);
    dhds[3]=0.25*(1-rvalue);


/*
for(i=0;i<2;i++)
{
  for(j=0;j<2;j++)
  printf("dhdr[%d][%d] %lf",i,j,dhdr[2*i+j] );
}printf("\n" );
*/
/*
for(i=0;i<2;i++)
{
  for(j=0;j<2;j++)
  printf("dhds[%d][%d] %lf",i,j,dhds[2*i+j] );
}printf("\n" );

for(i=0;i<2;i++)
{
  for(j=0;j<2;j++)
  printf("dhdr[%d][%d] %lf",i,j,dhdr[2*i+j] );
}printf("\n" );*/
}
void jacobi(double * jacobian, int nnel,int n, double *dhdr, double*dhds,double*xcoord, double *ycoord )
{
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++)
            jacobian[n*i+j]=0;
    }
    for(i=0;i<nnel;i++)
    {//printf("dhdr %lf   xc %lf\n",dhdr[i],xcoord[i] );
    jacobian[0]+=(dhdr[i]*xcoord[i]);//printf("J0= %lf= dhdr[%d]=%lf * xcord[%d]=%lf\n",jacobian[0],i, dhdr[i],i,xcoord[i] );
    jacobian[1]+=(dhdr[i]*ycoord[i]);//printf("J1= %lf= dhdr[%d]=%lf * ycord[%d]=%lf\n",jacobian[1],i, dhdr[i],i,ycoord[i] );
    jacobian[2]+=(dhds[i]*xcoord[i]);//printf("J2= %lf= dhds[%d]=%lf * xcord[%d]=%lf\n",jacobian[2],i, dhds[i],i,xcoord[i] );
    jacobian[3]+=(dhds[i]*ycoord[i]);//printf("J3= %lf= dhds[%d]=%lf * ycord[%d]=%lf\n",jacobian[3],i, dhds[i],i,ycoord[i] );
    }
}
/*
double det (double *kk, int n)
{
  int i,j,k; double ratio;  double a[n][n];

  for(i = 0; i < n; i++)
       for(j = 0; j < n; j++)
       a[i][j]=kk[n*i+j];

  for(i = 0; i < n; i++){
       for(j = 0; j < n; j++){
           if(j>i){
               ratio = a[j][i]/a[i][i];
               for(k = 0; k < n; k++){
                   a[j][k] -= ratio * a[i][k];
               }
           }
       }
     }
     double detj;
     detj=1;
     for(i=0;i<n;i++)
     detj *= a[i][i];
     return detj;
   }
*/
double det(double *a, int n)
{
    double detj;int i,j;
    /*  for(i=0;i<n;i++)
    {  for(j=0;j<n;j++)
    printf("a %lf",a[n*i+j] );
    printf("\n" );
    }
    */
    detj=(a[0]*a[3])-(a[1]*a[2]);
    //printf("detj %lf\n",detj );
    return detj;
}
void invert(double*a, int n,double det, double *ai)
{

    double adj[n][n];int i,j;
    adj[0][0]=a[3];
    adj[0][1]=-a[1];
    adj[1][0]=-a[2];
    adj[1][1]=a[0];


    for(i=0;i<n;i++)
    { 
        for(j=0;j<n;j++)
        {
            ai[n*i+j]=(adj[i][j]/det);
        }
    }
//printf("det  %lf\n",det );
/*
for(i=0;i<n;i++)
{  for(j=0;j<n;j++)
printf("a[%d][%d] %lf",i,j,a[n*i+j] );
printf("\n" );
}printf("\n" );

for(i=0;i<n;i++)
{  for(j=0;j<n;j++)
printf("ai[%d][%d] %lf",i,j,ai[n*i+j] );
printf("\n" );
}printf("\n" );
*/
}

/*void invert(double *a,int n, double det, double *aa);
{
  double adjoint[n][n];int i,j,k; double cof;
  if(det==0)
  {return;}
for(i=0;i<n;i++)
{
  for(j=0;j<n;j++)
  {
    for(k=0;k<n;k++)
  }
}
}
*/

void derivshape(double *dhdx,double *dhdy,int n,int nnel,double *dhdr,double *dhds,double *jinv)
{
    int i;
    for(i=0;i<nnel;i++)
    {
        //printf("jinv  %lf\n",jinv[1] );
        dhdx[i]=jinv[n*0+0]*dhdr[i]+jinv[n*0+1]*dhds[i];
        dhdy[i]=jinv[n*1+0]*dhdr[i]+jinv[n*1+1]*dhds[i];
    }
}

void bmat(double *be,int nnel,int n, double *dhdx,double *dhdy)
{
    // here not parametrized!!
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<nnel;j++)
        {
            if(i==0)
            {
                be[nnel*i+j]=dhdx[j];
            }
            if(i==1)
            {
                be[nnel*i+j]=dhdy[j];
            }
        }
    }
}

void btdb(double*kb, double *be, double wtx, double wty, double detj,int n, int nnel)
{
    int i,j,k; double res[nnel][nnel]; double bet[nnel][n];
    for(i=0;i<nnel;i++)
    {
        for(j=0;j<n;j++)
        {
            bet[i][j]=be[nnel*j+i];
        }
    }
    //  printf("detj %lf\n",detj );
    for(i=0;i<nnel;i++)
    {
        for(j=0;j<nnel;j++)
            res[i][j]=0;
    }
    for(i=0;i<nnel;i++)
    {
        for(j=0;j<nnel;j++)
        {
            for(k=0;k<n;k++)
            {
                res[i][j]+=(bet[i][k]*be[nnel*k+j]);
            }//printf("res[%d][%d]  = %lf",i,j,res[i][j] );
        }//printf("\n" );
    }//printf("\n" );

    /*  for(i=0;i<nnel;i++)
    {  for(j=0;j<nnel;j++)
        printf("kbefore[%d][%d] %lf",i,j,kb[nnel*i+j] );
        printf("\n" );
    }printf("\n" );
    */
    for(i=0;i<nnel;i++)
    {
        for(j=0;j<nnel;j++)
        {
            kb[nnel*i+j]+=((res[i][j])*wtx*wty*detj);
        //  printf("kb[%d][%d] =  %lf",i,j,kb[nnel*i+j] );
        }//printf("\n" );
    }
}



void quad2(int nglx,int ngly,int mnode,int ngl, double * gqp, double *gqw )
{
    double px[nglx];double py[ngly];
    double wx[nglx];double wy[ngly];int i;
    for(i=0;i<nglx;i++)
        wx[i]=0;px[i]=0;
    for(i=0;i<ngly;i++)
        wy[i]=0;py[i]=0;
    quad1(px,wx,nglx);
    quad1(py,wy,ngly);
    for(i=0;i<nglx;i++)
    {
        gqp[mnode*i+0]=px[i];
        gqw[mnode*i+0]=wx[i];
    }
    for(i=0;i<ngly;i++)
    {
        gqp[mnode*i+1]=py[i];
        gqw[mnode*i+1]=wy[i];
    }

}

//stiffmatgen(nnel,sysdof,iel,nel,nd,p,K);

void stiffmatgen(int n,int N, int iel, int nel, int *nd, double* k, double * K )
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            K[N*(nd[i]-1)+(nd[j]-1)]+=k[n*i+j];
        // printf("K[%d] [%d] = %lf  k %lf ",(nd[i]-1),(nd[j]-1),  K[N*(nd[i]-1)+(nd[j]-1)],k[n*i+j]);
        }//printf("\n" );
    }
}

void applycond (double *k, double *f,int n,int m, int * bcdof, double * bcval )
{
    int i,j,a;
    for(i=0;i<n;i++)
    {
        a=bcdof[i]-1;
        for(j=0;j<m;j++)
        {
            k[a*m+j]=0;
        }
        k[a*m+a]=1;
        f[a]=bcval[i];
    }
}

int main()
{
    char cwd[512] = __FILE__;
    int null_index;
    // To find the parent directory of the source file
    for(int i = 511;i >= 0;i--){
        if(cwd[i] == '\0'){
            null_index = i;
            break;
        }
    }
    #ifdef _WIN32
    for(int j = null_index-1;j >= 0;j--){
        if(cwd[j] == '\\'){
            break;
        }
        else{
            cwd[j] = '\0';
        }
    }
    #elif __linux__
    for(int j = null_index-1;j >= 0;j--){
        if(cwd[j] == '/'){
            break;
        }
        else{
            cwd[j] = '\0';
        }
    }
    #endif
    int i,j,k,nnode,mnode,nel,mel;
    double *coord;
    int *elem;
    k=0;i=0;j=0;
    FILE *f1;FILE *f3;
    FILE *f2;FILE *f4;
    char nodeinfo[512];
    char elementinfo[512];
    strcpy(nodeinfo,cwd);
    strcat(nodeinfo,"nodeinfo.txt");
    strcpy(elementinfo,cwd);
    strcat(elementinfo,"eleminfo.txt");
    f1=fopen(nodeinfo,"r");
    if(f1==NULL)
    {
        printf("failed to open nodeinfo matrix\n" );
        exit;
    }
    else
        printf("succesfully opened nodeinfo matrix\n" );
    while(!feof(f1)&&k<2)
    {
        if(k==0)
        {
            fscanf(f1, "%d",&nnode);
            k++;
        }
        if(k==1)
        {
            fscanf(f1,"%d",&mnode);
            k++;
        }
    }
    fclose(f1);

    k=0;
    f2=fopen(elementinfo,"r");
    if(f2==NULL)
    {
        printf("failed to open eleminfo matrix\n" );
        exit;
    }
    else
        printf("succesfully opened eleminfo matrix\n" );
    while(!feof(f2)&&k<2)
    {
        if(k==0)
        {
            fscanf(f2, "%d",&nel);
            k++;
        }
        if(k==1)
        {
            fscanf(f2,"%d",&mel);
            k++;
        }
    }
    k=0;fclose(f2);
    coord=(double*)malloc(nnode*mnode*sizeof(double));
    elem=(int*)malloc(nel*mel*sizeof(int));
    char prefix[255];
    char elements[255];
    char nodes[255];
    printf("Enter the prefix of the files:\n");
    scanf("%s",&prefix);
    strcpy(elements,cwd);
    strcpy(nodes,cwd);
    strcat(elements,prefix);
    strcat(nodes,prefix);
    strcat(elements,"_elements.txt");
    strcat(nodes,"_nodes.txt");
    printf("Elements file: %s, Nodes file: %s\n",elements,nodes);
    f3=fopen(nodes,"r");
    if(f3==NULL)
    {
        printf("failed to open unstrucquad_nodes \n" );
        exit;
    }
    else
        printf("succesfully opened unstrucquad_nodes \n" );
    while(!feof(f3))
    {
        for(i=0;i<nnode;i++)
        {
            for(j=0;j<mnode;j++)
            {
                fscanf(f3, "%lf",&coord[mnode*i+j]);
            }
        }
    }
    fclose(f3);
    for(i=0;i<nnode;i++)
    {
        for(j=0;j<mnode;j++)
            printf("%lf ",coord[mnode*i+j] );
        printf("\n" );
    }
    f4=fopen(elements,"r");
    if(f4==NULL)
    {
        printf("failed to open unstrucquad_elements \n" );
        exit;
    }
    else
        printf("succesfully opened unstrucquad_elements \n" );
    while(!feof(f4))
    {
        for(i=0;i<nel;i++)
        {
            for(j=0;j<mel;j++)
            {
                fscanf(f4, "%d",&elem[mel*i+j]);
            }
        }
    }
    fclose(f4);
    for(i=0;i<nel;i++)
    {
        for(j=0;j<mel;j++)
            printf("%d ",elem[mel*i+j] );
        printf("\n" );
    }
    int nnel=4;int m,l;
    int ngqx,ngqy; int ngl;
    ngqx=2;ngqy=2;
    if(ngqx>ngqy)
        ngl=ngqx;
    else
        ngl=ngqy;
    int dofpernode=1;int nbc;int s;
    int tst=0;
    for(i=0;i<nnode;i++)
    {
        if(coord[mnode*i+0]==0&&coord[mnode*i+1]!=0&&coord[mnode*i+1]!=1)
        {
            tst++;
        }
        if(coord[mnode*i+0]==1&&coord[mnode*i+1]!=0&&coord[mnode*i+1]!=1)
        {
            tst++;
        }
        if(coord[mnode*i+1]==0)
        {
            tst++;
        }
        if(coord[mnode*i+1]==1)
        {
            tst++;
        }
    }
    //printf("test %d",tst);
    nbc=tst;
    double*gqw; double*gqp;double*ycoord;double*xcoord; int quadp; int totdof;
    quadp=4;
    int *nd; double *be;double *shape; double *dhdr; double * dhds;
    double *jacobian; double * jinv; double *temp;double *dhdx;double *dhdy;
    totdof=nnel*1; double *p; double *K;int *bcnode; double *bcval;
    double *F;int iel;double qx,qy,wtx,wty,detj;
    int sysdof=nnode*dofpernode;
    gqw=(double*)malloc(ngl*mnode*sizeof(double));
    gqp=(double*)malloc(ngl*mnode*sizeof(double));
    nd=(int*)malloc(nnel*sizeof(int));
    xcoord=(double*)malloc(nnel*sizeof(double));
    ycoord=(double*)malloc(nnel*sizeof(double));
    shape=(double*)malloc(nnel*sizeof(double));
    dhdr=(double*)malloc(nnel*sizeof(double));
    dhds=(double*)malloc(nnel*sizeof(double));
    jacobian=(double*)malloc(mnode*mnode*sizeof(double));
    jinv=(double*)malloc(mnode*mnode*sizeof(double));
    temp=(double*)malloc(mnode*mnode*sizeof(double));
    dhdx=(double*)malloc(nnel*sizeof(double));
    dhdy=(double*)malloc(nnel*sizeof(double));
    //NOT PARAMETRIZED!!
    be=(double*)malloc(4*2*sizeof(double));
    p=(double*)malloc(nnel*nnel*sizeof(double));
    K=(double*)malloc(sysdof*sysdof*sizeof(double));
    bcnode=(int*)malloc(nbc*sizeof(int));
    bcval=(double*)malloc(nbc*sizeof(double));
    F=(double*)malloc(sysdof*sizeof(double));

    // -----------------------------------------------------------------------------
    // **** IMPORTANT SECTION : Boundary Conditions ****
    for(i=0;i<nbc;i++)
    {
        bcnode[k]=i+1;
        bcval[k]=0;
    }
    k=0;
    for(i=0;i<nnode;i++)
    {
        if(coord[mnode*i+0]==0 && coord[mnode*i+1]!=0 && coord[mnode*i+1]!=1) // u(0,y) for all 0 < y < 1
        {
            bcnode[k]=i+1;
            bcval[k]=0;
            k++;
        }
        if(coord[mnode*i+0]==1 && coord[mnode*i+1]!=0 && coord[mnode*i+1]!=1) // u(1,y) for all 0 < y < 1
        {
            bcnode[k]=i+1;
            bcval[k]=0;
            k++;
        }
        if(coord[mnode*i+1]==0) // u(x,0) for all 0 <= x <= 1
        {
            bcnode[k]=i+1;
            bcval[k]=0;
            k++;
        }
        if(coord[mnode*i+1]==1) // u(x,1) for all 0 <= x <= 1
        {
            bcnode[k]=i+1;
            bcval[k]=1;
            k++;
        }
    }
    // -----------------------------------------------------------------------------
    //printf("\n k  %d, nbc   %d",k,nbc );
    k=0;i=0;iel=0;
    for(i=0;i<sysdof;i++)
    {
        F[i]=0;
        for(j=0;j<sysdof;j++)
        {
            K[sysdof*i+j]=0;
            //printf("sys  %d\n",sysdof*i+j );
        }
    }
    for(i=0;i<nnel;i++)
        for(j=0;j<mnode;j++)
            be[mnode*i+j]=0;
    quad2(ngqx,ngqy,mnode,ngl,gqp,gqw);

    /*
    for(i=0;i<ngqx;i++)
    {
    for(j=0;j<mnode;j++)
    //  gqp[mnode*i+j]=gqp[mnode*j+i];
    //printf("point[%d][%d] %lf   weight[%d][%d] %lf\n",i,j,gqp[mnode*i+j],i,j,gqw[mnode*i+j] );
    }
    */

    for(iel=0;iel<nel;iel++)
    {
        for(i=0;i<nnel;i++)
        {//printf("iel = %d\n",iel );

            //printf("i  %d\n",i );
            nd[i]=elem[mel*iel+i];
            //printf("nd  %d\n",nd[i] );
            xcoord[i]=coord[(mnode*(nd[i]-1))+0];
            ycoord[i]=coord[(mnode*(nd[i]-1))+1];
            //  printf("xcoord = %lf  ycoord %lf \n",xcoord[i],ycoord[i] );
        }
        for(i=0;i<nnel;i++)
            for(j=0;j<nnel;j++)
                p[nnel*i+j]=0;
        for(i=0;i<ngqx;i++)
        {
            qx=gqp[mnode*i+0];
            wtx=gqw[mnode*i+0];
            for(j=0;j<ngqy;j++)
            {
                qy=gqp[mnode*j+1];
                wty=gqw[mnode*j+1];
                //  printf("wy %lf qy %lf\n",wty,qy );
                shapefn(qx,qy,shape,dhdr,dhds);
                /*  for(k=0;k<mnode;k++)
                {for(m=0;m<mnode;m++)
                    printf("shape %lf",shape[mnode*k+m] );
                    printf("\n" );
                }
                */
                jacobi(jacobian,nnel,mnode,dhdr,dhds,xcoord,ycoord);
                /*
                for(k=0;k<mnode;k++)
                {for(m=0;m<mnode;m++)
                    printf("J %lf",jacobian[mnode*k+m] );
                    printf("\n" );
                }*/

                /*
                for(i=0;i<4;i++)
                printf("dhds[%d] %lf\n",i,dhds[i] );
                for(i=0;i<4;i++)
                printf("dhdr[%d] %lf\n",i,dhdr[i] );

                */
                detj=det(jacobian,mnode);
                //printf("%lf\n",detj );
                invert(jacobian,mnode,detj,jinv);
                /*for(k=0;k<mnode;k++)
                {
                for(m=0;m<mnode;m++)
                printf("%lf",jinv[mnode*i+j] );
                printf("\n" );
                }
                */
                derivshape(dhdx,dhdy,mnode,nnel,dhdr,dhds,jinv);
                bmat(be,nnel,mnode,dhdx,dhdy);
                /*
                for(k=0;k<mnode;k++)
                {for(m=0;m<nnel;m++)
                printf("be %lf ",be[nnel*k+m] );
                printf("\n" );
                }
                */

                btdb(p,be,wtx,wty,detj,mnode,nnel);
                /*
                for(k=0;k<nnel;k++)
                {for(m=0;m<nnel;m++)
                printf(" %lf ",p[mnode*k+m] );
                printf("\n" );
                }  printf("\n" );
                */
            }
        }
        stiffmatgen(nnel,sysdof,iel,nel,nd,p,K);

    }
    applycond(K,F,nbc,sysdof,bcnode,bcval);
    printf(" \n n = %d\n",sysdof );
    /*for(k=0;k<sysdof;k++)
    {for(m=0;m<sysdof;m++)
    {
    printf("%lf ",K[sysdof*k+m] );
    }
    printf("    %lf\n",F[k] );
    printf("\n" );
    }
    */
    FILE *o1,*o2,*o3,*o4;
    char kinfo[512];
    char kmat[512];
    char fvec[512];
    strcpy(kinfo,cwd);
    strcat(kinfo,"kinfo.txt");
    strcpy(kmat,cwd);
    strcat(kmat,"kmat.txt");
    strcpy(fvec,cwd);
    strcat(fvec,"Fvec.txt");
    o1=fopen(kinfo,"w");
    fprintf(o1, "%d\n",sysdof );
    fclose(o1);
    o2=fopen(kmat,"w");
    for(k=0;k<sysdof;k++)
    {
        for(m=0;m<sysdof;m++)
            fprintf(o2, "%lf\n",K[sysdof*k+m] );
    }
    fclose(o2);
    o3=fopen(fvec,"w");
    for(k=0;k<sysdof;k++)
    fprintf(o3, "%lf\n",F[k] );
    fclose(o3);
    free(coord);
    free(xcoord);
    free(ycoord);
    free(nd);
    free(elem);
    free(gqp);
    free(gqw);
    free(dhds);
    free(dhdr);
    free(shape);
    free(jacobian);
    free(jinv);
    free(temp);
    free(dhdx);
    free(dhdy);
    free(K);
    free(F);
    free(bcnode);
    free(bcval);
    free(p);
    free(be);
    return 0;
}
