/***
 * millipede: FastMarching.h
 * Added by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_FASTMARCHING
#define H_MILLIPEDE_FASTMARCHING

#include <cfloat>
#include <cmath>
#include <list>

#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_integral.hpp>

#include <itkCastImageFilter.h>
#include <itkConstantBoundaryCondition.h>
#include <itkConstShapedNeighborhoodIterator.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <common/adts/PriorityQueue.h>
#include <common/util/ITKImageUtil.h>
#include <common/util/NullType.h>
#include <common/util/QuadEqn.h>

namespace mp {

template <unsigned int Dimension, typename InputPixelType = unsigned int>
class FastMarching
{
	//#################### TEMPLATE PARAMETER CONSTRAINTS ####################
	BOOST_MPL_ASSERT_MSG(boost::is_integral<InputPixelType>::value,
						 NON_INTEGRAL_PIXEL_TYPES_ARE_NOT_SUPPORTED,
						 (InputPixelType));

	//#################### TYPEDEFS ####################
private:
	typedef double TimePixelType;	

	typedef itk::Image<short, Dimension> GradientMagnitudeImage;
	typedef itk::Image<InputPixelType, Dimension> InputImage;
	typedef itk::Image<TimePixelType, Dimension> TimeImage;

	typedef typename InputImage::IndexType Index;
	typedef std::list<Index> Indices;

	typedef typename GradientMagnitudeImage::Pointer GradientMagnitudeImagePointer;
	typedef typename InputImage::Pointer InputImagePointer;
	typedef typename TimeImage::Pointer TimeImagePointer;

	typedef PriorityQueue<Index, TimePixelType, NullType, std::less<TimePixelType>, itk::Functor::IndexLexicographicCompare<Dimension> > PQ;

	typedef itk::Offset<Dimension> NeighbourOffset;
	typedef std::vector<NeighbourOffset> NeighbourOffsets;

	typedef itk::ConstShapedNeighborhoodIterator<TimeImage> ConstShapedNeighbourhoodIterator;
	typedef typename ConstShapedNeighbourhoodIterator::ConstIterator ConstNeighbourIterator;
	typedef itk::ConstNeighborhoodIterator<TimeImage> ConstNeighbourhoodIterator;
		
	//#################### PRIVATE VARIABLES ####################
private:
	GradientMagnitudeImagePointer m_gradients;
	TimeImagePointer m_times;
	NeighbourOffsets m_offsets;

	//#################### CONSTRUCTORS ####################
public:
	FastMarching(const InputImagePointer& input, const Indices& initial)
	{
		construct_gradients(input);
		run(initial);
	}

	FastMarching(const GradientMagnitudeImagePointer& gradients, const Indices& initial)
	:	m_gradients(gradients)	
	{
		run(initial);
	}	

	//#################### PUBLIC METHODS ####################
public:
	Indices get_shape_at_time(TimePixelType time)
	{
		Indices indices;
		itk::ImageRegionConstIterator<TimeImage> it(m_times, m_times->GetLargestPossibleRegion());
		for(it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			if(it.Get() <= time)
			{
				indices.push_back(it.GetIndex());
			}
		}

		return indices;
	}

	Indices get_shape_at_first_stop(int threshold)
	{
		typedef std::list<TimePixelType> Values;
		// Sort time image values in increasing order by inserting all of them into a priority queue and then popping out. 
		Indices allIndices;
		Values allValues;

		PQ pq;
		itk::ImageRegionConstIteratorWithIndex<TimeImage> it(m_times, m_times->GetLargestPossibleRegion());
		for(it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			pq.insert(it.GetIndex(), it.Get(), NullType());
		}
		allIndices.resize(pq.size());
		allValues.resize(pq.size());
		while(!pq.empty())
		{
			typename PQ::Element e = pq.top();
			pq.pop();
			allIndices.push_back(e.id());
			allValues.push_back(e.key());
		}

		// Repeatedly increase the time value starting from 0 by a fixed delta until the expansion of the region is sufficiently small. For more information refer to my dissertation.
		TimePixelType delta = 2.0;
		typename Indices::iterator curIndexIt = allIndices.begin();
		Values::iterator curValueIt = allValues.begin();
		Values::iterator largestValueIt = allValues.end();
		int difference = 0;
		for(TimePixelType value=delta; ; value+=delta)
		{
			difference = 0;
			while(curValueIt != largestValueIt && *curValueIt <= value)
			{
				++curIndexIt;
				++curValueIt;
				++difference;
			}

			if(difference < threshold || curValueIt == largestValueIt)	break;
		}

		// Take the indices with currentRegionSize-many smallest values (this is the region for time value) and return it.
		Indices indices;
		indices.splice(indices.begin(), allIndices, allIndices.begin(), curIndexIt);
		return indices;
	}

	//#################### PRIVATE METHODS ####################
private:
	static PQ build_initial_propagation_queue(const Indices& initial)
	{
		PQ pq;
		for(typename Indices::const_iterator jt=initial.begin(), jend=initial.end(); jt!=jend; ++jt)
		{
			pq.insert(*jt, 0.0, NullType());	
		}	

		return pq;
	}

	TimePixelType compute_new_time_at(const Index& index, const ConstNeighbourhoodIterator& it)
	{		
		double A, B, C;

		double alpha = 1.0;
		double coefficient = 10.0;	
		A = 0;
		B = 0;
		// The following formulas are used: C = -1 / F^2 , where F(x) = coefficient * e^(-2*alpha*x). For more information refer to my dissertation.
		C = -(exp(2*alpha*m_gradients->GetPixel(index))) / (coefficient*coefficient);
		TimePixelType time = it.GetCenterPixel();
		for(unsigned int i=0; i<Dimension; ++i)
		{
			// The offsets in the i-th direction have indices Dimension-1-i and Dimension+i in the vector of offsets.
			TimePixelType adjTime = std::min(it.GetPixel(m_offsets[Dimension-1-i]), it.GetPixel(m_offsets[Dimension+i]));
			if(time > adjTime)
			{
				// The spacing in the time image is considered to be 1. 
				A += 1;
				B += -2 * adjTime;
				C += adjTime * adjTime;
			}	
		}

		return QuadEqn::largest_root(A, B, C);
	}	

	void construct_gradients(const InputImagePointer& input)
	{
		typedef itk::Image<float, Dimension> RealImage;
		typedef typename RealImage::Pointer RealImagePointer;

		// Cast the input image to make its pixels real-valued.
		typedef itk::CastImageFilter<InputImage,RealImage> CastFilter;
		typename CastFilter::Pointer castFilter = CastFilter::New();
		castFilter->SetInput(input);
		castFilter->Update();
		RealImagePointer real = castFilter->GetOutput();

		// Smooth this real image using anisotropic diffusion filtering.
		typedef itk::GradientAnisotropicDiffusionImageFilter<RealImage,RealImage> ADFilter;
		typename ADFilter::Pointer adFilter = ADFilter::New();
		adFilter->SetInput(real);
		adFilter->SetConductanceParameter(1.0);
		adFilter->SetNumberOfIterations(20);
		adFilter->SetTimeStep(0.0625);
		adFilter->Update();
		real = adFilter->GetOutput();

		// Calculate the gradient magnitude of the smoothed image.
		typedef itk::GradientMagnitudeImageFilter<RealImage, GradientMagnitudeImage> GMFilter;
		typename GMFilter::Pointer gmFilter = GMFilter::New();
		gmFilter->SetInput(real);
		gmFilter->SetUseImageSpacingOff();
		gmFilter->Update();
		m_gradients = gmFilter->GetOutput();
	}	

	void initialise_times()
	{
		m_times = TimeImage::New();
		m_times->SetRegions(m_gradients->GetLargestPossibleRegion());
		m_times->Allocate();

		itk::ImageRegionIterator<TimeImage> it(m_times, m_times->GetLargestPossibleRegion());
		for(it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			it.Set(DBL_MAX);
		}
	}

	void propagate_surface(PQ& pq)
	{
		TimePixelType largeNumber = 100.0;

		itk::Size<Dimension> radius;
		radius.Fill(1);
		ConstShapedNeighbourhoodIterator it(radius, m_times, m_times->GetLargestPossibleRegion());

		for(typename NeighbourOffsets::const_iterator jt=m_offsets.begin(), jend=m_offsets.end(); jt!=jend; ++jt)
		{
			it.ActivateOffset(*jt);	
		}

		itk::ConstantBoundaryCondition<TimeImage> itCondition;
		itCondition.SetConstant(0.0);
		it.OverrideBoundaryCondition(&itCondition);	

		ConstNeighbourhoodIterator kt(radius, m_times, m_times->GetLargestPossibleRegion());
		itk::ConstantBoundaryCondition<TimeImage> ktCondition;
		ktCondition.SetConstant(DBL_MAX);
		kt.OverrideBoundaryCondition(&ktCondition);

		while(!pq.empty())
		{
			typename PQ::Element e = pq.top();
			pq.pop();
			Index cur = e.id();
			TimePixelType curTime = e.key();

			// To avoid real-valued type precision errors stop the process when it reaches some large time value.
			if(curTime > largeNumber)	break;

			m_times->SetPixel(cur, curTime);
			it.SetLocation(cur);
			for(ConstNeighbourIterator jt=it.Begin(), jend=it.End(); jt!=jend; ++jt)
			{
				Index adj = cur + jt.GetNeighborhoodOffset();
				TimePixelType adjTime = jt.Get();

				if(adjTime <= curTime)	continue;

				kt.SetLocation(adj);
				TimePixelType adjNewTime = compute_new_time_at(adj, kt);
				if(pq.contains(adj))
				{
					if (adjNewTime < pq.element(adj).key())
					{
						pq.update_key(adj, adjNewTime);
					}
				}
				else
				{
					pq.insert(adj, adjNewTime, NullType());
				}
			}
		}		
	}
	
	void run(const Indices& initial)
	{
		m_offsets = ITKImageUtil::make_2d_connected_offsets<Dimension>();

		initialise_times();

		PQ pq =  build_initial_propagation_queue(initial);
		propagate_surface(pq);
	}	
};

}

#endif
