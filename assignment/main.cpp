#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle, vtkRenderingFreeType, vtkRenderingOpenGL2)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2)

#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageData.h>
#include <vtkImageActor.h>
#include <vtkImageViewer2.h>
#include <vtkCamera.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <math.h>
#include "Vector3.h"


/*
	Reads the file and store the data in a vector
	@param filename The name of the file
	@param width The width of the data (obtained from the first line in the file)
	@param height The height of the data (obtained from the first line in the file)
	@return A vector with the data. Each cell is a Vector3 pointer
*/
std::vector<Vector3*>* readDataFile(const char* filename, int &width, int &height) 
{
	// Initialize the file reader
	std::ifstream infile(filename);

	// Initialize the vector where data will be stored
	std::vector<Vector3*>* data = new std::vector<Vector3*>();

	// Read the first line and get the width and height of the data
	std::string line;
	std::getline(infile, line);
	std::istringstream iss(line);
	double a, b, c;
	iss >> a >> b >> c;
	width = (int)a;
	height = (int)b;

	// Read the file lines
	while (std::getline(infile, line))
	{
		// Read the line, generate the new vector3 and store it in the vector
		iss.clear();
		iss.str(line);
		iss >> a >> b >> c;
		data->push_back(new Vector3(a, b, c));
	}

	// Close the file and return the data
	infile.close();
	return data;
}


/*
	Generates a width x height random noise image
	@param width The width (in pixels) of the image
	@param height The height (in pixels) of the image
	@return The pointer to the image with the random noise
*/
vtkSmartPointer<vtkImageData> generateWhiteNoise(int width, int height)
{
	// Initialize the imgae where the noise will be stored
	vtkSmartPointer<vtkImageData> noise = vtkSmartPointer<vtkImageData>::New();
	noise->SetDimensions(width, height, 1);
	noise->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	// Initialize the random number generator and the distribution where the random numbers will come from
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, 255);

	// Traverse the image by pixel
	for (int y = 0; y < height; y += 1)
	{
		for (int x = 0; x < width; x += 1)
		{
			// Get an uniform random variable
			int randomValue = distribution(generator);

			// Get the pointer to the current pixel and set the random value
			unsigned char* pixel = static_cast<unsigned char*>(noise->GetScalarPointer(x, y, 0));
			pixel[0] = (unsigned char)(randomValue);
			pixel[1] = (unsigned char)(randomValue);
			pixel[2] = (unsigned char)(randomValue);
		}
	}

	// Return the pointer to the noise image
	return noise;
}


/*
	Calculates a new point using bilinear interpolation
	@param P1 The bottom left point
	@param P2 The bottom right point
	@param P3 The top right point
	@param P4 The top left point
	@param u The horizontal parameter
	@param v The vertical parameter
	@return The interpolated vector
*/
Vector3* bilinearInterpolation(Vector3* P1, Vector3* P2, Vector3* P3, Vector3* P4, double u, double v)
{
	// Calculate the location of the point
	Vector3* term1 = P1->clone();
	Vector3* term2 = P2->clone()->sub(P1)->multiply(u);
	Vector3* term3 = P4->clone()->sub(P1)->multiply(v);
	Vector3* term4 = P1->clone()->sub(P2)->add(P3)->sub(P4)->multiply(u * v);
	Vector3* result = term1->clone()->add(term2)->add(term3)->add(term4);

	// Delete pointers
	delete term1;
	delete term2;
	delete term3;
	delete term4;

	// Return the result
	return result;
}


/*
	Calculates the interpolated 2D vector for a (x, y) pixel using the data
	@param x The x coordinate of the pixel
	@param y The y coordinate of the pixel
	@param data The vector with the values of the data grid
	@param n The width of the data grid
	@param m The height of the data grid
	@return A Vector3 pointer with the interpolated value. This is actually a 2D vector with z = 0
*/
Vector3* getPixel2DVector(int x, int y, std::vector<Vector3*>* data, int n, int m)
{
	// Calculate the (u, v) values of the pixel in terms of the
	double u = ((double)(x * (n - 1))) / 299.0;
	double v = ((double)(y * (m - 1))) / 299.0;

	// Get the coordinates for the four closest points in the data grid
	int u_min = floor(u);
	int v_min = floor(v);
	int u_max = ceil(u);
	int v_max = ceil(v);

	// Get the four points for the bilinear interpolation
	Vector3* P1 = data->at(u_min + (v_min * n));
	Vector3* P2 = data->at(u_max + (v_min * n));
	Vector3* P3 = data->at(u_max + (v_max * n));
	Vector3* P4 = data->at(u_min + (v_max * n));

	// Get (u, v) for the bilinear interpolation
	// NOTE: They must be values between [0, 1]
	u = u - u_min;
	v = v - v_min;

	// Calculate the new vector using bilinear interpolation and return it
	return bilinearInterpolation(P1, P2, P3, P4, u, v);
}


/*
	Fourth order Runge-Kutta (RK4) for numerical integration
	@param x The x coordinate of the starting point
	@param y The y coordinate of the starting point
	@param h The step size for the next coordinate in the streamline
	@param data The vector with the values of the data grid
	@param n The width of the data grid
	@param m The height of the data grid
	@return A Vector3 pointer with the coordinates of the next pixel in the streamline
*/
Vector3* RungeKutta4(int x, int y, int width, int height, double h, std::vector<Vector3*>* data, int n, int m)
{
	// Define the vector with the coordinates of the current point
	Vector3* Xn = new Vector3(x, y, 0);

	// Calculate the K components:
	Vector3* K1 = getPixel2DVector(x, y, data, n, m)->normalize();

	Vector3* K2 = NULL;
	int xn = x + (((double)h / 2.0) * K1->x);
	int yn = y + (((double)h / 2.0) * K1->y);
	if (xn >= 0 && xn < width && yn >= 0 && yn < height) 
	{
		K2 = getPixel2DVector(xn, yn, data, n, m)->normalize();
	}
	else 
	{
		return new Vector3(-1, -1, 0);
	}

	Vector3* K3 = NULL;
	xn = x + (((double)h / 2.0) * K2->x);
	yn = y + (((double)h / 2.0) * K2->y);
	if (xn >= 0 && xn < width && yn >= 0 && yn < height)
	{
		K3 = getPixel2DVector(xn, yn, data, n, m)->normalize();
	}
	else
	{
		return new Vector3(-1, -1, 0);
	}

	Vector3* K4 = NULL;
	xn = x + ((double)h * K3->x);
	yn = y + ((double)h * K3->y);
	if (xn >= 0 && xn < width && yn >= 0 && yn < height)
	{
		K4 = getPixel2DVector(xn, yn, data, n, m)->normalize();
	}
	else
	{
		return new Vector3(-1, -1, 0);
	}

	// Calculate the next coordinate
	Vector3* term1 = K1->clone();
	Vector3* term2 = K2->clone()->multiply(2.0);
	Vector3* term3 = K3->clone()->multiply(2.0);
	Vector3* term4 = K4->clone();
	Vector3* Xnext = term1->clone()->add(term2)->add(term3)->add(term4)->multiply(h / 6.0)->add(Xn);
	
	// Delete pointers
	delete Xn;
	delete K1;
	delete K2;
	delete K3;
	delete K4;
	delete term1;
	delete term2;
	delete term3;
	delete term4;

	// Return the result
	return Xnext;
}


/*
	Returns the noise value stored in the given (x, y) pixel in the given noise image
	@param x The x coordinate of the pixel
	@param y The y coordinate of the pixel
	@param noise The pointer to the noise image
	@return The noise value stored in the given pixel
*/
double getNoiseValue(int x, int y, vtkSmartPointer<vtkImageData> noise)
{
	// Get the pointer to the current pixel in the noise image
	unsigned char* pixel = static_cast<unsigned char*>(noise->GetScalarPointer(x, y, 0));

	// Get the noise value and return it
	// NOTE: The RGB channels have the same noise value. Then, just return the value stored in
	// the red channel
	int color = pixel[0];
	return color;
}


/*
	Returns the average value of the elements in the given vector
	@param values A pointer to the vector with the values
	@return The average value of the elements in the vector
*/
double average(std::vector<Vector3*>* values, vtkSmartPointer<vtkImageData> noise, int width, int height)
{
	// Get the number of elements in the vector
	size_t n = values->size();

	// Initialize the sum of the elements
	double sum = 0.0;
	int N = 0;

	// Traverse the elements of the vector and add them up
	for (int i = 0; i < n; i += 1) 
	{
		Vector3* p = values->at(i);
		int x = p->x;
		int y = p->y;

		if (x >= 0 && x < width && y >= 0 && y < height) 
		{
			double value = getNoiseValue(p->x, p->y, noise);
			sum += value;

			N += 1;
		}
	}

	// Return the average value
	return (N > 0) ? sum / ((double)N) : 0.0;
}


/*
	Generates a linear gaussian kernel of length 2L - 1. L is use as the standard deviation value
	@param L The length of each tai of the gaussian (from the center to the end)
	@return A vector with the values of the gaussian
*/
std::vector<double>* getGaussianKernel(int L) 
{
	// Set the value of pi
	double pi = 3.14159265359;

	// Define the standard deviation value
	double E = (double)L / 2.0;

	// Intialize the vector where the gaussian values will be stored
	std::vector<double>* gaussian = new std::vector<double>();

	// Generate the gaussian values
	for (int l = -L + 1; l < L; l += 1) 
	{
		// Generate the gaussian value for the current value of l and insert it into the 
		// gaussian vector
		double value = (1.0 / (E * sqrt(2 * pi))) * exp(-((l * l) / (2.0 * E * E)));
		gaussian->push_back(value);
	}

	// Return the vector with the gaussian values
	return gaussian;
}


/*
	Returns a vector with the coordinates of the elements of the streamline that start in the given pixel
	@param x The x coordinate of the pixel
	@param y The y coordinate of the pixel
	@param width The width of the requested image
	@param height The height of the requested image
	@param L The length of each streamline section
	@param h The step length for the Runge-Kutta4 algorithm
	@param data The vector with the data grid information
	@param n The width of the data grid
	@param m The height of the data grid
	@param noise The pointer to the noise image
	@return The pointer to the vector with the streamline coordinates
*/
std::vector<Vector3*>* getStreamLine(int x, int y, int width, int height, int L, double h, std::vector<Vector3*>* data, int n, int m, vtkSmartPointer<vtkImageData> noise)
{
	// Initialize the vector where the noise values of the positive streamline will be stored
	std::vector<Vector3*>* positiveStreamline = new std::vector<Vector3*>();

	// Set the initial location of the positive streamline
	int x0 = x;
	int y0 = y;

	// Traverse through the positive section of the streamline
	for (int l = 0; l < L; l += 1) 
	{
		// Store the coordinates of the current pixel
		positiveStreamline->push_back(new Vector3(x0, y0, 0));

		// If the current location is a valid location then find the next one. Otherwise set (-1, -1)
		if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < width) 
		{
			// Get the location of the next positive streamline pixel
			Vector3* next = RungeKutta4(x0, y0, width, height, h, data, n, m);
			x0 = next->x;
			y0 = next->y;
			delete next;
		}
		else 
		{
			x0 = -1;
			y0 = -1;
		}
	}

	// Repeat the same steps for the negative streamline section

	// Initialize the vector where the noise values of the negative streamline will be stored
	std::vector<Vector3*>* negativeStreamline = new std::vector<Vector3*>();

	// Set the initial location of the negative streamline
	x0 = x;
	y0 = y;

	// Traverse through the negative section of the streamline
	for (int l = 0; l > -L; l -= 1)
	{
		// Store the coordinates of the current pixel
		negativeStreamline->push_back(new Vector3(x0, y0, 0));

		if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height) 
		{
			// Get the location of the next negative streamline pixel
			Vector3* next = RungeKutta4(x0, y0, width, height, -h, data, n, m);
			x0 = next->x;
			y0 = next->y;
			delete next;
		}
		else 
		{
			x0 = -1;
			y0 = -1;
		}
	}

	// Now we have both sections of the streamline that fall within the region of the requested image. Merge both 
	// sections of the streamline into a single vector

	// Initialize the vector where the information of the streamline will be stored
	std::vector<Vector3*>* streamline = new std::vector<Vector3*>();

	// Insert the elements of the negative streamline sections
	// NOTE: Insert them in the inverse order for having the values in the right order
	for (int i = L - 1; i >= 0; i -= 1) 
	{
		streamline->push_back(negativeStreamline->at(i));
	}

	// Insert the elements of the positive streamline
	// NOTE: Do not insert the first element since it was already inserted from the negative streamline
	for (int i = 1; i < L; i += 1) 
	{
		streamline->push_back(positiveStreamline->at(i));
	}

	// Delete the vector with the positive and negative streamline sections (we don't need them anymore)
	delete negativeStreamline;
	delete positiveStreamline;

	// Return the vector with the streamline values
	return streamline;
}


/*
	Linear Integral Convolution
	@param width The width of the requested image
	@param height The height of the requested image
	@param L The length of each streamline section
	@param h The step for traversing the streamlines
	@param data The vector with the values of the data grid
	@param n The width of the data grid
	@param m The height of the data grid
	@return A vtkSmartPointer<vtkImageData> pointer to the resultant image
*/
vtkSmartPointer<vtkImageData> LIC(int width, int height, int L, double h, std::vector<Vector3*>* data, int n, int m)
{
	// Initialize the image
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
	image->SetDimensions(width, height, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	// Get a white noise image with the same dimensions as the requested image
	vtkSmartPointer<vtkImageData> noise = generateWhiteNoise(width, height);

	// Get the gaussian kernel
	int length = (2 * L) - 1;
	std::vector<double>* gaussian = getGaussianKernel(L);

	// Traverse the image by pixel
	for (int y = 0; y < height; y += 1)
	{
		for (int x = 0; x < width; x += 1)
		{
			// Get the streamline associated with the current pixel
			std::vector<Vector3*>* streamline = getStreamLine(x, y, width, height, L, h, data, n, m, noise);

			// Get the average noise value found in the streamline
			double avg = average(streamline, noise, width, height);

			// Initialize the sum values
			double sum_num = 0.0;
			double sum_den = 0.0;

			// Traverse the elements of the streamline
			// NOTE: The length of the streamline should be 2L - 1
			for (int i = 0; i < length; i += 1) 
			{
				// Get the coordinates of the current pixel in the streamline
				Vector3* p = streamline->at(i);
				int x0 = p->x;
				int y0 = p->y;

				// If the coordinates fall within the region of the image then process it
				if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height) 
				{
					sum_num += (getNoiseValue(x0, y0, noise) * gaussian->at(i));
					sum_den += gaussian->at(i);
					//sum_num += (getNoiseValue(x0, y0, noise) * avg);
					//sum_den += avg;
					
				}
			}
			
			// Using the obtained information calculate the color for the current pixel
			double value = sum_num / sum_den;
			if (value > 255) 
			{
				std::cout << value << std::endl;
			}

			int color = (int)value;
			
			// Get the pointer to the current pixel and set the random value
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			pixel[0] = (unsigned char)(color);
			pixel[1] = (unsigned char)(color);
			pixel[2] = (unsigned char)(color);
		}
	}

	// Return the image
	return image;
}


/*
	The main function
*/
int main(int argc, char** argv)
{
	// Read the file selection and define the file name
	int option = 0;
	std::cout << "1 = circle, 2 = wind: ";
	std::cin >> option;
	const char* filename = (option == 1) ? "circle.txt" : "wind.txt";

	int L;
	int h;
	std::cout << "L = ";
	std::cin >> L;
	std::cout << "h = ";
	std::cin >> h;

	// Read the data file and keep its dimensions (n = width, m = height)
	int n, m;
	std::vector<Vector3*>* data = readDataFile(filename, n, m);
	
	// Generate the LIC for visualizing the vector field
	vtkSmartPointer<vtkImageData> image = LIC(300, 300, L, h, data, n, m);
	//vtkSmartPointer<vtkImageData> image = generateWhiteNoise(300, 300);

	// Define the renderer window interactor and start the program
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

	// Set the image viewer
	vtkSmartPointer<vtkImageViewer2> viewer = vtkSmartPointer<vtkImageViewer2>::New();
	viewer->SetInputData(image);
	viewer->SetupInteractor(renderWindowInteractor);
	viewer->GetRenderer()->ResetCamera();
	viewer->GetRenderer()->GetActiveCamera()->SetParallelScale(150);
	viewer->Render();

	// Setup the render window interactor and start the program
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}
