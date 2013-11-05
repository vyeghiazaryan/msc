/***
 * millipede: LiverIdentifier3D.cpp
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#include "LiverIdentifier3D.h"

#include <climits>

#include <boost/bind.hpp>

#include <common/dicom/volumes/DICOMVolume.h>
#include <common/interfacemotion/FastMarching.h>
#include <common/util/ITKImageUtil.h>

namespace mp {

//#################### CONSTRUCTORS ####################
LiverIdentifier3D::LiverIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF)
:	FeatureIdentifier(dicomVolume, volumeIPF)
{}

//#################### PUBLIC METHODS ####################
int LiverIdentifier3D::length() const
{
	return 6;
}

//#################### PRIVATE METHODS ####################
void LiverIdentifier3D::execute_impl()
{
	set_status("Identifying liver...");

	VolumeIPFMultiFeatureSelection_Ptr multiFeatureSelection = get_multi_feature_selection();

	// Step 1: Filter for the liver candidates.
	std::list<PFNodeID> candidates = filter_branch_nodes(boost::bind(&LiverIdentifier3D::is_candidate, this, _1, _2));
	increment_progress();

	// Step 2: Pick the best candidate (namely, the one which stretches furthest to the left of the image).
	PFNodeID bestCandidate;
	int bestXMin = INT_MAX;
	for(std::list<PFNodeID>::const_iterator it=candidates.begin(), iend=candidates.end(); it!=iend; ++it)
	{
		const BranchProperties& properties = volume_ipf()->branch_properties(*it);
		if(properties.x_min() < bestXMin)
		{
			bestCandidate = *it;
			bestXMin = properties.x_min();
		}
	}

	// If we can't find a candidate liver, exit.
	if(bestCandidate == PFNodeID::invalid()) return;

	increment_progress();	

	// Step 3: Convert the partition forest candidate nodes to positions in the volume.
	std::list<itk::Index<3> > positions;
	std::deque<int> leaves = volume_ipf()->receptive_region_of(bestCandidate);
	for(std::deque<int>::const_iterator jt=leaves.begin(), jend=leaves.end(); jt!=jend; ++jt)
	{
		if(160 < volume_ipf()->leaf_properties(*jt).grey_value() && volume_ipf()->leaf_properties(*jt).grey_value() < 180)
		{
			positions.push_back(volume_ipf()->position_of_leaf(*jt));
		}
	}	

	increment_progress();

	// Step 4: Use the Fast Marching Method.
	FastMarching<3> fm(volume_ipf()->leaf_layer()->gradient_magnitude_image(), positions);
	int slices = dicom_volume()->size()[2];
	positions = fm.get_shape_at_first_stop(3000 * slices);

	increment_progress();

	// Step 5: Convert the positions in the volume to partition forest selection.
	std::set<int> leafIndices;
	for(std::list<itk::Index<3> >::const_iterator it=positions.begin(), iend=positions.end(); it!=iend; ++it)
	{
		leafIndices.insert(volume_ipf()->leaf_of_position(*it));
	}
	PartitionForestSelection_Ptr filledRegion(new PartitionForestSelectionT(volume_ipf(), leafIndices));

	increment_progress();

	// Step 6: Mark the result as liver.
	multiFeatureSelection->identify_selection(filledRegion, AbdominalFeature::LIVER);
}

bool LiverIdentifier3D::is_candidate(const PFNodeID& node, const BranchProperties& properties) const
{
	itk::Index<3> volumeSize = ITKImageUtil::make_index_from_size(dicom_volume()->size());
	int sliceCount = properties.z_max() + 1 - properties.z_min();

	return	node.layer() == 1 &&								// it should be in the lowest branch layer (higher nodes generally stretch too far)
			properties.voxel_count() >= MIN_VOXELS_PER_SLICE * sliceCount &&	// it should be reasonably-sized
			MIN_MEAN_GREY_VALUE <= properties.mean_grey_value() && 
			properties.mean_grey_value() <= MAX_MEAN_GREY_VALUE &&			// it should have a reasonable grey value
			properties.x_min() < volumeSize[0]/2;					// and at least part of it should be on the left-hand side of the image
}

}
