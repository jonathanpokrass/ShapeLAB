#include <mex.h>

#include "pcl/common/common_headers.h"
#include "pcl/features/normal_3d.h"
#include "pcl/kdtree/kdtree_flann.h"
#include "pcl/features/fpfh.h"
#include "pcl/features/pfh.h"

#include "boost/numeric/ublas/vector.hpp"

using namespace boost::numeric::ublas;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgTxt
	 *       within an if statement, because it will never get to the else
	 *       statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
	 *       the MEX-file) 
	*/
	if(nrhs < 2){ 
		mexErrMsgTxt("At least two inputs required.");
	}
	if(nrhs > 3){
		mexErrMsgTxt("At most three inputs required.");
	}
	if(nlhs < 2){
		mexErrMsgTxt("Seven output required.");
	}

	/* check to make sure the first input argument is a square matrix with double entry */
	if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
		mexErrMsgTxt("Input x must be a real value matrix.");
	}
	mwSize np = mxGetN(prhs[0]);
	mwSize dim = mxGetM(prhs[0]);

	double *ctdim = mxGetPr(prhs[1]);
    if(mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1){
		mexErrMsgTxt("Second input must be a scalar.");
	}
	double radius = (ctdim[0]);

	double* points;
	points = mxGetPr(prhs[0]);

	
	// Create point cloud
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::Normal>::Ptr cloudNormals(new pcl::PointCloud<pcl::Normal>);
	double* curPoint = points;

	for(int p = 0; p < np; p++)
	{
		cloud->points.push_back(pcl::PointXYZ(curPoint[0], curPoint[1], curPoint[2]));
		//cloudNormals->points.push_back(pcl::Normal(1,0,0));
		curPoint += 3;
	}
	cloud->width = (int) cloud->points.size();
	cloud->height = 1; 
	
	// Estimate normals
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
	pcl::search::KdTree<pcl::PointXYZ>::Ptr cloudSearchTree(new pcl::search::KdTree<pcl::PointXYZ>());
	

	ne.setInputCloud(cloud);
	ne.setSearchMethod(cloudSearchTree);
	ne.setKSearch(10);
	//ne.setRadiusSearch(radius * 0.7);
	
	ne.compute(*cloudNormals); 
	
	// Calculate features
	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh;
	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfh_features(new pcl::PointCloud<pcl::FPFHSignature33>());
	fpfh.setSearchMethod(cloudSearchTree);
	fpfh.setInputCloud(cloud);
	fpfh.setRadiusSearch(radius);
	//fpfh.setKSearch(15);
	fpfh.setInputNormals(cloudNormals);
	fpfh.compute(*fpfh_features);

	//
	if(nlhs > 2)
	{

		pcl::PFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::PFHSignature125> pfh;
		pcl::PointCloud<pcl::PFHSignature125>::Ptr pfh_features(new pcl::PointCloud<pcl::PFHSignature125>());
		pfh.setSearchMethod(cloudSearchTree);
		pfh.setInputCloud(cloud);
		pfh.setRadiusSearch(radius);
		//fpfh.setKSearch(15);
		pfh.setInputNormals(cloudNormals);
		pfh.compute(*pfh_features);
		plhs[2] = mxCreateDoubleMatrix(125, np, mxREAL);
		
		double* featureData = mxGetPr(plhs[2]);

		for (int i = 0; i < np; i++)
		{
			const int FEATURE_SIZE = 125;
			for (int j = 0; j < FEATURE_SIZE; j++)
			{
				featureData[j] = pfh_features->points[i].histogram[j];
			}
			featureData += FEATURE_SIZE;
		}
	}


	// Return results

	mxArray* normArr = NULL;
	double* normData = NULL;
	double* featureData = NULL;
	if(nlhs > 1)
	{
		normArr = mxCreateDoubleMatrix(3, np, mxREAL); 
		plhs[1] = normArr;
		normData = mxGetPr(normArr);
	} 

	plhs[0] = mxCreateDoubleMatrix(33, np, mxREAL);
	featureData = mxGetPr(plhs[0]);

	for (int i = 0; i < np; i++)
	{
		normData[0] = cloudNormals->points[i].normal_x;
		normData[1] = cloudNormals->points[i].normal_y;
		normData[2] = cloudNormals->points[i].normal_z;
		normData+=3;
		const int FEATURE_SIZE = 33;
		for (int j = 0; j < FEATURE_SIZE; j++)
		{
			featureData[j] = fpfh_features->points[i].histogram[j];
		}
		featureData += FEATURE_SIZE;
	}


	if(nlhs > 2)
	{
	
	}
	
}


