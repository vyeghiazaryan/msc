/***
 * millipede: DICOMImageLeafLayer.cpp
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#include "DICOMImageLeafLayer.h"

namespace mp {

//#################### CONSTRUCTORS ####################
DICOMImageLeafLayer::DICOMImageLeafLayer(const std::vector<DICOMPixelProperties>& nodeProperties, int sizeX, int sizeY, int sizeZ)
{
	initialise(nodeProperties, sizeX, sizeY, sizeZ);
}

DICOMImageLeafLayer::DICOMImageLeafLayer(const itk::Image<int,2>::Pointer& baseImage,
										 const itk::Image<unsigned char,2>::Pointer& windowedImage,
										 const itk::Image<short,2>::Pointer& gradientMagnitudeImage)
{
	assert(baseImage->GetLargestPossibleRegion().GetSize() == windowedImage->GetLargestPossibleRegion().GetSize());

	typedef itk::Image<int,2> BaseImage;
	typedef itk::Image<short,2> GradientMagnitudeImage;
	typedef itk::Image<unsigned char,2> WindowedImage;

	const BaseImage::SizeType& size = baseImage->GetLargestPossibleRegion().GetSize();

	std::vector<DICOMPixelProperties> nodeProperties;
	nodeProperties.reserve(size[0] * size[1]);

	BaseImage::IndexType baseIndex;
	GradientMagnitudeImage::IndexType gradientMagnitudeIndex;
	WindowedImage::IndexType windowedIndex;
	for(size_t y=0; y<size[1]; ++y)
		for(size_t x=0; x<size[0]; ++x)
		{
			baseIndex[0] = gradientMagnitudeIndex[0] = windowedIndex[0] = x;
			baseIndex[1] = gradientMagnitudeIndex[1] = windowedIndex[1] = y;
			nodeProperties.push_back(DICOMPixelProperties(baseImage->GetPixel(baseIndex),
														  gradientMagnitudeImage->GetPixel(gradientMagnitudeIndex),
														  windowedImage->GetPixel(windowedIndex)));
		}

	initialise(nodeProperties, size[0], size[1]);
}

DICOMImageLeafLayer::DICOMImageLeafLayer(const itk::Image<int,3>::Pointer& baseImage,
										 const itk::Image<unsigned char,3>::Pointer& windowedImage,
										 const itk::Image<short,3>::Pointer& gradientMagnitudeImage)
{
	assert(baseImage->GetLargestPossibleRegion().GetSize() == windowedImage->GetLargestPossibleRegion().GetSize());

	typedef itk::Image<int,3> BaseImage;
	typedef itk::Image<short,3> GradientMagnitudeImage;
	typedef itk::Image<unsigned char,3> WindowedImage;

	const BaseImage::SizeType& size = baseImage->GetLargestPossibleRegion().GetSize();

	std::vector<DICOMPixelProperties> nodeProperties;
	nodeProperties.reserve(size[0] * size[1] * size[2]);

	BaseImage::IndexType baseIndex;
	GradientMagnitudeImage::IndexType gradientMagnitudeIndex;
	WindowedImage::IndexType windowedIndex;
	for(size_t z=0; z<size[2]; ++z)
		for(size_t y=0; y<size[1]; ++y)
			for(size_t x=0; x<size[0]; ++x)
			{
				baseIndex[0] = gradientMagnitudeIndex[0] = windowedIndex[0] = x;
				baseIndex[1] = gradientMagnitudeIndex[1] = windowedIndex[1] = y;
				baseIndex[2] = gradientMagnitudeIndex[2] = windowedIndex[2] = z;
				nodeProperties.push_back(DICOMPixelProperties(baseImage->GetPixel(baseIndex),
															  gradientMagnitudeImage->GetPixel(gradientMagnitudeIndex),
															  windowedImage->GetPixel(windowedIndex)));
			}

	initialise(nodeProperties, size[0], size[1], size[2]);
}

//#################### PUBLIC METHODS ####################
// Precondition: has_edge(u, v)
DICOMImageLeafLayer::EdgeWeight DICOMImageLeafLayer::edge_weight(int u, int v) const
{
	return std::max(m_nodes[u].properties().gradient_magnitude_value(), m_nodes[v].properties().gradient_magnitude_value());
}

itk::Image<short,3>::Pointer DICOMImageLeafLayer::gradient_magnitude_image() const
{
	typedef itk::Image<short,3> Image;
	
	Image::Pointer image = Image::New();
	Image::SizeType size;
	size[0] = m_sizeX;
	size[1] = m_sizeY;
	size[2] = m_sizeZ;
	
	Image::RegionType region;
	region.SetSize(size);
	image->SetRegions(region);
	image->Allocate();

	Image::IndexType index;
	for(size_t z=0; z<size[2]; ++z)
		for(size_t y=0; y<size[1]; ++y)
			for(size_t x=0; x<size[0]; ++x)
			{
				index[0] = x;
				index[1] = y;
				index[2] = z;
				int n = x + y * m_sizeX + z * m_sizeXY;
				image->SetPixel(index, m_nodes[n].properties().gradient_magnitude_value());
			}

	return image;
}

}
