#include "mex.h"
#include <math.h>

// % calculate the depth of an non-integer point via bilinear interpolation
//  v1 ------- v2
//   |  .p     |
//   |         |
//   v4-------v3
double depth_Bilinear(double px, double py, int imgm, int imgn, double *depth_img)
{
    int zx, zy, cidx;
    double coeff1, coeff2, coeff3, coeff4;
    double depth1, depth2, depth3, depth4;
    double depth_pts = 0;
    zx = (int)floor(px);
    zy = (int)floor(py);

    if (zx>0 && zx<imgn && zy>0 && zy<imgm)
    {
        coeff1 = (zx+1-px)*(zy+1-py);
        coeff2 = (px-zx)*(zy+1-py);
        coeff3 = (px-zx)*(py-zy);
        coeff4 = (zx+1-px)*(py-zy);
        cidx = (zx-1)*imgm + zy-1; 
        depth1 = depth_img[cidx];
        depth2 = depth_img[cidx+imgm];
        depth3 = depth_img[cidx+imgm+1];
        depth4 = depth_img[cidx+1];
        depth_pts = coeff1*depth1 + coeff2*depth2 + coeff3*depth3 + coeff4*depth4;
    }

    return depth_pts;
}


/* give a point and a quadrangle, find whether it is inside the quadrangle*/
bool JudgeInsideQuad( double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4, double px, double py)
{
    double Denominator1=0, Denominator2=0;
    double Numerator1=0, Numerator2=0, Numerator3=0;
    double Numerator4=0, Numerator5=0, Numerator6=0;
    bool InsideFlag1 = false;
    bool InsideFlag2 = false;

    Denominator1 = (px2-px1)*(py3-py1)-(px3-px1)*(py2-py1);
    Numerator1   = (px2-px)*(py3-py)-(px3-px)*(py2-py);
    Numerator2   = (px-px1)*(py3-py1)-(px3-px1)*(py-py1);
    Numerator3   = (px2-px1)*(py-py1)-(px-px1)*(py2-py1);
    
    if ((Denominator1>0 && Numerator1>0 && Numerator2>0 && Numerator3>0) || (Denominator1<0 && Numerator1<0 && Numerator2<0 && Numerator3<0))
    {
        InsideFlag1 = true; 
    }
    else
    {
        Denominator2 = (px3-px1)*(py4-py1)-(px4-px1)*(py3-py1);
        Numerator4   = (px3-px)*(py4-py)-(px4-px)*(py3-py);
        Numerator5   = (px-px1)*(py4-py1)-(px4-px1)*(py-py1);
        Numerator6   = (px3-px1)*(py-py1)-(px-px1)*(py3-py1);
        
        if ((Denominator2>0 && Numerator4>0 && Numerator5>0 && Numerator6>0) || (Denominator2<0 && Numerator4<0 && Numerator5<0 && Numerator6<0))
        {
            InsideFlag2 = true; 
        }
    }

    return InsideFlag1 || InsideFlag2; 

}

/* give a point and a quadrangle, find the bilinear interpolation coefficients alpha, beta*/
/*% Cal_AlphaInQuad 此函数计算四边形内部的一个点如何被四边形的四个顶点双线性表示出来，假设此点在四边形中（内部包括边界）
% 具体公式：P=(1-alpha)*(1-beta)*P1+alpha*(1-beta)*P2+alpha*beta*P3+(1-alpha)*beta*P4
% （四边形顶点顺序为P1与P3相对，P2与P1纵坐标接近，处在相近的水平线上）*/
void Cal_AlphaBetaInQuad(double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4, double px, double py, double bi_coeff[2])
{	
	/* version-2  by liaotianli using maple 18*/ 
	double Alpha=0, Beta=0;
	double BetaCoeff_a=0, BetaCoeff_b=0, BetaCoeff_c=0;
	double Betadelta = 0, Beta1=0, Beta2=0;
	
	BetaCoeff_a = (px1-px4)*(py2-py3)-(py1-py4)*(px2-px3);
	BetaCoeff_b = (-py1+py2-py3+py4)*px+(px1-px2+px3-px4)*py-2*px1*py2+px1*py3+2*px2*py1-px2*py4-px3*py1+px4*py2;
	BetaCoeff_c = (py1-py2)*px+(-px1+px2)*py+px1*py2-px2*py1;
		
	if (fabs(BetaCoeff_a)<=1e-15)
	{	
		if (fabs(BetaCoeff_b)>1e-15)  //  to avoid system crush
		{
			Beta = -BetaCoeff_c/BetaCoeff_b;
		}
	}
	else
	{
        Betadelta = BetaCoeff_b*BetaCoeff_b-4*BetaCoeff_a*BetaCoeff_c;
        if (Betadelta>=0)  //  to avoid system crush
        {
            Beta1 = (-BetaCoeff_b + sqrt(Betadelta))/(2*BetaCoeff_a);
            Beta2 = (-BetaCoeff_b - sqrt(Betadelta))/(2*BetaCoeff_a);
            if (Beta1>=0 && Beta1<=1)
            {
                Beta = Beta1;
            }
            else
            {
                Beta = Beta2;
            }
        }
	}
    
    if (fabs(Beta*px1-Beta*px2+Beta*px3-Beta*px4-px1+px2)>1e-15)  //  to avoid system crush
    {
        Alpha = (Beta*px1-Beta*px4+px-px1)/(Beta*px1-Beta*px2+Beta*px3-Beta*px4-px1+px2);
    }
	   
    if (Alpha>=0 && Alpha<=1 && Beta>=0 && Beta<=1)
    {
        bi_coeff[0] = Alpha;
        bi_coeff[1] = Beta;
    }
    else
    {
        bi_coeff[0]=0; bi_coeff[1]=0;  // to avoid system crush
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Input/output variables. */
    double *img1; /* original input image. */
    double *img2; /* original input image. */
    double *depth_img1;  /* depth map of img1*/
    double *depth_img2;  /* depth map of img2*/
    double dch;  /*  height of warped image */
    double dcw;  /*  width of warped image */
    double *X;      /* X positions of the original mesh grid. */
    double *Y;      /* Y positions of the original mesh grid. */ 
    double *wX;     /* X positions of the warped mesh grid. */
    double *wY;     /* Y positions of the warped mesh grid. */   
    double *off;   /* offset of warped image r.t. original image */
    double *epipolar_He; /* infinity homography and epipole of img2->img1*/
    double *occlude_grid; /* occlusion index in the non-overlapping region */
    double para_dist; /* parameter of displacement*/
              
    double *warped_img1;     /* Warped img (img after being warped with mesh transformation */
    double *warped_mask1;
    double *warped_errors;
    double *warped_depth1;

    /* Intermediate variables.*/
    int ch, cw, C1, C2;  /* int-type ch, cw, C1, C2*/
    int imgm, imgn;      /* 'm' and 'n': size of the original image. */    
    int x, y;            /* x,y positions in each warped grid*/
    int mleft, mright, mup, mdown;   /* the external rectangle of each warped grid */

	int indij;          /* 1D-indice of X,Y (or wX, wY)*/
    int i, j;          

    bool InsideFlag;        /* judge whether a point is inside quadrangle */
	double alpha, beta;   /* alpha, beta coefficients of bilinear representation */
	double bi_coeff[2]={0};
    double posx, posy;
    int px, py, cidx, sidx;        /*contain the indexes of the warped and original images respectively */
    double tmpdepth, map_x, map_y;
    double dist_disp;
    double tmperror, tmpr, tmpg, tmpb;
    double tmp_grid;

    /* Since MATLAB works with 2D matrices and C works with 1D arrays we need the following variables for generating the 
     * final stitched RGB 2D matrix for MATLAB, this image/matrix contains the warped img2. */
    int ch1canv; /* Displacement for channel 1 in image 1 (R). */
    int ch2canv; /* Displacement for channel 2 in image 1 (G). */
    int ch3canv; /* Displacement for channel 3 in image 1 (B). */
    
    int ch1img; /* Displacement for channel 1 in image 2 (R). */
    int ch2img; /* Displacement for channel 2 in image 2 (G). */
    int ch3img; /* Displacement for channel 3 in image 2 (B). */
    
    
//     /* Check for proper number of arguments. */    
//     if (nrhs!=15)
//     {
//         mexErrMsgTxt("Wrong number of inputs.");
//     }
//     else if (nlhs!=1)
//     {
//         mexErrMsgTxt("Wrong number of output arguments.");
//     }
// 	if (!mxIsDouble(prhs[0]))
// 	{
// 		mexErrMsgTxt("Input must be type of double.");
// 	}
	    
    /* Assign pointers to inputs. */   
    img1=mxGetPr(prhs[0]);
    img2=mxGetPr(prhs[1]);
    depth_img1=mxGetPr(prhs[2]);
    depth_img2=mxGetPr(prhs[3]);
    dch=mxGetScalar(prhs[4]);
    dcw=mxGetScalar(prhs[5]);
    X=mxGetPr(prhs[6]);
    Y=mxGetPr(prhs[7]);
    wX=mxGetPr(prhs[8]);
    wY=mxGetPr(prhs[9]); 
    off=mxGetPr(prhs[10]);  
    epipolar_He=mxGetPr(prhs[11]);
    occlude_grid=mxGetPr(prhs[12]);
    para_dist=mxGetScalar(prhs[13]);

    imgm=mxGetM(prhs[0]);
    imgn=mxGetN(prhs[0])/3; /* It is an RGB image, which means that for MATLAB it is an Mx(Nx3) image. So, we have  
                                  * to divide N by 3 in order to get the proper size of the image in C. */
    C1=mxGetM(prhs[6])-1; /* rows of mesh grid */
    C2=mxGetN(prhs[6])-1; /* columns of mesh grid */
    ch=(int)dch;  /* height of warped image */
    cw=(int)dcw;  /* width of warped image */                        
    
    /* Create matrix for the return arguments. */
    plhs[0]=mxCreateDoubleMatrix( ch, cw*3, mxREAL);
    plhs[1]=mxCreateDoubleMatrix( ch, cw, mxREAL);
    plhs[2]=mxCreateDoubleMatrix( ch, cw, mxREAL);
    plhs[3]=mxCreateDoubleMatrix( ch, cw, mxREAL);

    /* Assign pointers to output canvas (warped image2). */
    warped_img1=mxGetPr(plhs[0]);
    warped_mask1=mxGetPr(plhs[1]);
    warped_errors=mxGetPr(plhs[2]);
    warped_depth1=mxGetPr(plhs[3]);       
    
    /* Initialize displacements. */
    ch1canv=0;
    ch2canv=ch*cw;
    ch3canv=ch*cw*2;
    
    ch1img=0;
    ch2img=imgn*imgm;
    ch3img=imgn*imgm*2;
    
    /*initialization*/
    for (j=0; j<cw;j++)
    {
        for (i=0;i<ch;i++)
        {
            warped_img1[j*ch+i+ch1canv]=0;
            warped_img1[j*ch+i+ch2canv]=0;
            warped_img1[j*ch+i+ch3canv]=0;
            warped_mask1[j*ch+i]=0;
            warped_errors[j*ch+i]=0;
            warped_depth1[j*ch+i]=0;
        }
    }
	
    
    /* Start computations. */
    /* For each grid. */
    for(j=0;j<C2;j++)
    {
        for(i=0;i<C1;i++)
        {            
			indij = j*(C1+1)+i;
            mleft=(int)floor(fmin(fmin(wX[indij], wX[indij+1]), fmin(wX[indij+C1+1], wX[indij+C1+2])));
            mright=(int)ceil(fmax(fmax(wX[indij], wX[indij+1]), fmax(wX[indij+C1+1], wX[indij+C1+2])));
            mup= (int)floor(fmin(fmin(wY[indij], wY[indij+1]), fmin(wY[indij+C1+1], wY[indij+C1+2])));
            mdown=(int)ceil(fmax(fmax(wY[indij], wY[indij+1]), fmax(wY[indij+C1+1], wY[indij+C1+2])));     
            /* Get grid pixel for current grid. */
            for(x=mleft; x<=mright; x++)
            {
                for(y=mup; y<=mdown; y++)
                {
                    /* use bilinear method to interpolate  */
                    InsideFlag=JudgeInsideQuad( wX[indij], wY[indij], wX[indij+C1+1], wY[indij+C1+1], wX[indij+C1+2], wY[indij+C1+2], wX[indij+1], wY[indij+1], x, y );
                    if (InsideFlag)
                    {
						Cal_AlphaBetaInQuad(wX[indij], wY[indij], wX[indij+C1+1], wY[indij+C1+1], wX[indij+C1+2], wY[indij+C1+2], wX[indij+1], wY[indij+1], x, y, bi_coeff);
						alpha=bi_coeff[0];
						beta=bi_coeff[1];
						posx=(1-alpha)*(1-beta)*X[indij] + alpha*(1-beta)*X[indij+C1+1] + alpha*beta*X[indij+C1+2] + (1-alpha)*beta*X[indij+1];
                        posy=(1-alpha)*(1-beta)*Y[indij] + alpha*(1-beta)*Y[indij+C1+1] + alpha*beta*Y[indij+C1+2] + (1-alpha)*beta*Y[indij+1];
                        px=(int)floor(posx);
                        py=(int)floor(posy);
                        cidx=(int)((x+off[0]-2)*ch+y+off[1]-2);
                        warped_mask1[cidx]=warped_mask1[cidx]+1;
                        if (x>=1 && x<=imgn && y>=1 && y<=imgm) // for pixels in the overlapping region
                        {
                            tmpdepth=depth_img2[(x-1)*imgm+y-1];
                            map_x=(epipolar_He[0]*x+epipolar_He[1]*y+epipolar_He[2]+epipolar_He[9]/tmpdepth)/(epipolar_He[6]*x+epipolar_He[7]*y+epipolar_He[8]+epipolar_He[11]/tmpdepth);
                            map_y=(epipolar_He[3]*x+epipolar_He[4]*y+epipolar_He[5]+epipolar_He[10]/tmpdepth)/(epipolar_He[6]*x+epipolar_He[7]*y+epipolar_He[8]+epipolar_He[11]/tmpdepth);
                            dist_disp=(map_x-posx)*(map_x-posx)+(map_y-posy)*(map_y-posy);
                            if (dist_disp>=para_dist*para_dist)
                            {
                                //mexPrintf("dist_disp=%f\n",dist_disp);
                                continue;
                            }
                            else if (px>=1 && px<imgn && py>=1 && py<imgm)
                            {
                                sidx=(int)(px-1)*imgm+py-1;
                                tmpr=(px+1-posx)*(py+1-posy)*img1[sidx+ch1img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch1img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch1img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch1img];
                                tmpg=(px+1-posx)*(py+1-posy)*img1[sidx+ch2img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch2img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch2img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch2img];
                                tmpb=(px+1-posx)*(py+1-posy)*img1[sidx+ch3img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch3img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch3img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch3img];
                                tmperror=(tmpr-img2[(x-1)*imgm+y-1])*(tmpr-img2[(x-1)*imgm+y-1])+(tmpg-img2[(x-1)*imgm+y-1+ch2img])*(tmpg-img2[(x-1)*imgm+y-1+ch2img])+(tmpb-img2[(x-1)*imgm+y-1+ch3img])*(tmpb-img2[(x-1)*imgm+y-1+ch3img]);
                                if (warped_mask1[cidx]<=1)
                                {
                                    warped_img1[cidx]=tmpr;
                                    warped_img1[cidx+ch2canv]=tmpg;
                                    warped_img1[cidx+ch3canv]=tmpb;
                                    warped_errors[cidx]=tmperror;
                                }
                                else if (tmperror<warped_errors[cidx])
                                {
                                    warped_img1[cidx]=tmpr;
                                    warped_img1[cidx+ch2canv]=tmpg;
                                    warped_img1[cidx+ch3canv]=tmpb;
                                    warped_errors[cidx]=tmperror;
                                }
                            }
                        }
                        else if (px>=1 && px<imgn && py>=1 && py<imgm)  /* if the current pixel in the target image falls inside canvas. */
                        {
                            tmp_grid=occlude_grid[j*C1+i];
                            if (tmp_grid==1)
                            {
                                continue;
                            }
                            tmpdepth=depth_Bilinear(posx, posy, imgm, imgn, depth_img1);
                            sidx=(int)(px-1)*imgm+py-1; 
                            if (warped_mask1[cidx]<=1)
                            {
                                tmpr=(px+1-posx)*(py+1-posy)*img1[sidx+ch1img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch1img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch1img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch1img];
                                tmpg=(px+1-posx)*(py+1-posy)*img1[sidx+ch2img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch2img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch2img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch2img];
                                tmpb=(px+1-posx)*(py+1-posy)*img1[sidx+ch3img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch3img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch3img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch3img];
                                warped_img1[cidx]=tmpr;
                                warped_img1[cidx+ch2canv]=tmpg;
                                warped_img1[cidx+ch3canv]=tmpb;
                                warped_depth1[cidx]=tmpdepth;
                            }
                            else if (tmpdepth<warped_depth1[cidx])
                            {
                                tmpr=(px+1-posx)*(py+1-posy)*img1[sidx+ch1img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch1img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch1img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch1img];
                                tmpg=(px+1-posx)*(py+1-posy)*img1[sidx+ch2img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch2img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch2img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch2img];
                                tmpb=(px+1-posx)*(py+1-posy)*img1[sidx+ch3img] + (px+1-posx)*(posy-py)*img1[sidx+1+ch3img] + (posx-px)*(py+1-posy)*img1[sidx+imgm+ch3img] + (posx-px)*(posy-py)*img1[sidx+imgm+1+ch3img];
                                warped_img1[cidx]=tmpr;
                                warped_img1[cidx+ch2canv]=tmpg;
                                warped_img1[cidx+ch3canv]=tmpb;
                                warped_depth1[cidx]=tmpdepth;
                            }                        
                        }
                    }
                }
            }
        }
    }
            
    /* End.*/
    return;
}
