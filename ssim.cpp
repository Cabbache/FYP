#include <iostream>
#include <opencv2/video.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>

using namespace std;
using namespace cv;

//From https://docs.opencv.org/4.x/dd/d3d/tutorial_gpu_basics_similarity.html
Scalar getMSSIM( const Mat& i1, const Mat& i2)
{
    const double C1 = 6.5025, C2 = 58.5225;
    /***************************** INITS **********************************/
    int d     = CV_32F;
    Mat I1, I2;
    i1.convertTo(I1, d);           // cannot calculate on one byte large values
    i2.convertTo(I2, d);
    Mat I2_2   = I2.mul(I2);        // I2^2
    Mat I1_2   = I1.mul(I1);        // I1^2
    Mat I1_I2  = I1.mul(I2);        // I1 * I2
    /*************************** END INITS **********************************/
    Mat mu1, mu2;   // PRELIMINARY COMPUTING
    GaussianBlur(I1, mu1, Size(11, 11), 1.5);
    GaussianBlur(I2, mu2, Size(11, 11), 1.5);
    Mat mu1_2   =   mu1.mul(mu1);
    Mat mu2_2   =   mu2.mul(mu2);
    Mat mu1_mu2 =   mu1.mul(mu2);
    Mat sigma1_2, sigma2_2, sigma12;
    GaussianBlur(I1_2, sigma1_2, Size(11, 11), 1.5);
    sigma1_2 -= mu1_2;
    GaussianBlur(I2_2, sigma2_2, Size(11, 11), 1.5);
    sigma2_2 -= mu2_2;
    GaussianBlur(I1_I2, sigma12, Size(11, 11), 1.5);
    sigma12 -= mu1_mu2;
    Mat t1, t2, t3;
    t1 = 2 * mu1_mu2 + C1;
    t2 = 2 * sigma12 + C2;
    t3 = t1.mul(t2);              // t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))
    t1 = mu1_2 + mu2_2 + C1;
    t2 = sigma1_2 + sigma2_2 + C2;
    t1 = t1.mul(t2);               // t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2))
    Mat ssim_map;
    divide(t3, t1, ssim_map);      // ssim_map =  t3./t1;
    Scalar mssim = mean( ssim_map ); // mssim = average of ssim map
    return mssim;
}

int main(int argc, char **argv){
	if (argc != 3){
		cerr << "Usage: " << argv[0] << " [video file 1] [video file 2]" << endl;
		return 1;
	}

	VideoCapture cap1(argv[1]);
	VideoCapture cap2(argv[2]);

	if(!cap1.isOpened()){
		cerr << "failed to open " << argv[1];
		return 2;
	}

	if(!cap2.isOpened()){
		cerr << "failed to open " << argv[2];
		return 2;
	}

	double average = 0;
	int count = 0;
	while (1){
		Mat frame1;
		Mat frame2;
		cap1 >> frame1;
		cap2 >> frame2;
		
		if (frame1.empty() != frame2.empty()){
			cerr << "Videos have different amount of frames" << endl;
			return 3;
		}

		if (frame1.empty()) break;
		
		Scalar mssim = getMSSIM(frame1, frame2);
		cerr << mssim << endl;
		average += mssim[0] + mssim[1] + mssim[2];
		count += 3;
	}

	cerr << "Average ssim: " << endl;
	cout << (average / count) << endl;
	return 0;
}
