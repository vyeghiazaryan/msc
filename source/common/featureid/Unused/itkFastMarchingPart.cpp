

#include <itkImage.h>
#include <itkFastMarchingImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <cmath>



{




typedef itk::Image<float, 3> RealImage;
typedef typename RealImage::Pointer RealImagePointer;

typedef itk::CastImageFilter<itk::Image<unsigned char, 3>,RealImage > CastFilter;
typename CastFilter::Pointer castFilter = CastFilter::New();
castFilter->SetInput(dicom_volume()->windowed_image(WindowSettings(40, 400)));
castFilter->Update();
RealImagePointer real = castFilter->GetOutput();

typedef itk::GradientAnisotropicDiffusionImageFilter<RealImage,RealImage> ADFilter;
for(int i=0; i<10; ++i)
{
	typename ADFilter::Pointer adFilter = ADFilter::New();
	adFilter->SetInput(real);
	adFilter->SetConductanceParameter(1.0);
	adFilter->SetNumberOfIterations(1);
	adFilter->SetTimeStep(0.0625);
	adFilter->Update();
	real = adFilter->GetOutput();
}

typedef itk::GradientMagnitudeImageFilter<RealImage, itk::Image<float, 3> > GMFilter;
typename GMFilter::Pointer gmFilter = GMFilter::New();
gmFilter->SetInput(real);
gmFilter->SetUseImageSpacingOff();
gmFilter->Update();
itk::Image<float, 3>::Pointer m_gradients = gmFilter->GetOutput();std::cerr<<"Gradients Computed";

itk::Image<double, 3>::Pointer speed;
double alpha = -1.0;
speed = itk::Image<double, 3>::New();
		speed->SetRegions(m_gradients->GetLargestPossibleRegion());
		speed->Allocate();

		itk::ImageRegionIterator<itk::Image<double, 3> > it(speed, speed->GetLargestPossibleRegion());
		itk::ImageRegionIterator<itk::Image<float, 3> > jt(m_gradients, m_gradients->GetLargestPossibleRegion());
		for(it.GoToBegin(); !it.IsAtEnd(); ++it, ++jt)
		{
			it.Set(exp(alpha * jt.Get()));
		}


 typedef itk::BinaryThresholdImageFilter< itk::Image<double, 3>, 
                        itk::Image<unsigned char, 3>   >    ThresholdingFilterType;
  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();

  thresholder->SetLowerThreshold( 0.0 );
  thresholder->SetUpperThreshold( 4000.0 );

  thresholder->SetOutsideValue(  0  );
  thresholder->SetInsideValue(  255 );

 typedef  itk::FastMarchingImageFilter< itk::Image<double, 3> ,
                              itk::Image<double, 3 > >    FastMarchingFilterType; 
 FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();


fastMarching->SetInput( speed );
  thresholder->SetInput( fastMarching->GetOutput() );


 typedef FastMarchingFilterType::NodeContainer           NodeContainer;
  typedef FastMarchingFilterType::NodeType                NodeType;
  NodeContainer::Pointer seeds = NodeContainer::New();
  //  Software Guide : EndCodeSnippet 
  

  itk::Image<unsigned char, 3>::IndexType  seedPosition;
  
  seedPosition[0] = 135;
  seedPosition[1] = 340;
  seedPosition[2] = 0;

  NodeType node;
  const double seedValue = 0.0;
  
  node.SetValue( seedValue );
  node.SetIndex( seedPosition );
  seeds->Initialize();
  seeds->InsertElement( 0, node );
  fastMarching->SetTrialPoints(  seeds  );


fastMarching->SetOutputSize( 
           dicom_volume()->windowed_image(WindowSettings(40, 400))->GetLargestPossibleRegion().GetSize());
  fastMarching->SetStoppingValue(  5000.0);

thresholder->Update();

for(unsigned int n = 0; n < 1/*NUMBER*/; ++n)
	{
		increment_progress();
		PartitionForestSelection_Ptr region(new PartitionForestSelectionT(volume_ipf()));
	
		itk::ImageRegionIteratorWithIndex<itk::Image<unsigned char, 3> > it(thresholder->GetOutput(), thresholder->GetOutput()->GetLargestPossibleRegion());
		for(it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if(it.Get() > 100)
			{
				int leafIndex = ((it.GetIndex())[2]*volumeSize[1]+(it.GetIndex())[1])*volumeSize[0]+(it.GetIndex())[0];
			region->select_node(PFNodeID(0, leafIndex));
			}
		}
		multiFeatureSelection->identify_selection(region, AbdominalFeature::Enum(AbdominalFeature::LEVEL1+n));
	}
}

