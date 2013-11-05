/***
 * millipede: SpleenIdentifier3D.cpp
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#include "SpleenIdentifier3D.h"

#include <climits>

#include <boost/bind.hpp>

#include <common/dicom/volumes/DICOMVolume.h>
#include <common/interfacemotion/FastMarching.h>
#include <common/util/ITKImageUtil.h>

namespace mp {

//#################### CONSTRUCTORS ####################
SpleenIdentifier3D::SpleenIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF)
:	FeatureIdentifier(dicomVolume, volumeIPF)
{}

//#################### PUBLIC METHODS ####################
int SpleenIdentifier3D::length() const
{
	return 3;
}

//#################### PRIVATE METHODS ####################
void SpleenIdentifier3D::execute_impl()
{
	set_status("Identifying spleen...");

	VolumeIPFMultiFeatureSelection_Ptr multiFeatureSelection = get_multi_feature_selection();

	// Step 1: Calculate the combined properties of all the nodes marked as part of the spine.
	BranchProperties spineProperties = multiFeatureSelection->properties_of(AbdominalFeature::VERTEBRA);

	// Step 2: Filter for spleen candidates.
	std::list<PFNodeID> candidates = filter_branch_nodes(boost::bind(&SpleenIdentifier3D::is_candidate, this, _1, _2, spineProperties));

	// Step 3: Pick the best candidate (namely, the one which stretches furthest to the right of the image).
	PFNodeID bestCandidate;
	int bestXMax = INT_MIN;
	for(std::list<PFNodeID>::const_iterator it=candidates.begin(), iend=candidates.end(); it!=iend; ++it)
	{
		const BranchProperties& properties = volume_ipf()->branch_properties(*it);
		if(properties.x_max() > bestXMax)
		{
			bestCandidate = *it;
			bestXMax = properties.x_max();
		}
	}

	// If we can't find a candidate spleen, exit.
	if(bestCandidate == PFNodeID::invalid()) return;

	// Step 4: Convert the partition forest candidate nodes to positions in the volume.
	std::list<itk::Index<3> > positions;
	std::deque<int> leaves = volume_ipf()->receptive_region_of(bestCandidate);
	for(std::deque<int>::const_iterator jt=leaves.begin(), jend=leaves.end(); jt!=jend; ++jt)
	{
		if(140 < volume_ipf()->leaf_properties(*jt).grey_value() && volume_ipf()->leaf_properties(*jt).grey_value() < 170)
		{
			positions.push_back(volume_ipf()->position_of_leaf(*jt));
		}
	}

	increment_progress();

	// Step 5: Use the Fast Marching Method.
	FastMarching<3> fm(volume_ipf()->leaf_layer()->gradient_magnitude_image(), positions);
	positions = fm.get_shape_at_time(5.0);

	increment_progress();

	// Step 6: Convert the positions in the volume to partition forest selection.
	std::set<int> leafIndices;
	for(std::list<itk::Index<3> >::const_iterator it=positions.begin(), iend=positions.end(); it!=iend; ++it)
	{
		leafIndices.insert(volume_ipf()->leaf_of_position(*it));
	}
	PartitionForestSelection_Ptr region(new PartitionForestSelectionT(volume_ipf(), leafIndices));

	// Step 7: Mark the results as spleen.
	multiFeatureSelection->identify_selection(region, AbdominalFeature::SPLEEN);
}

bool SpleenIdentifier3D::is_candidate(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties) const
{
	itk::Index<3> volumeSize = ITKImageUtil::make_index_from_size(dicom_volume()->size());
	int sliceCount = properties.z_max() + 1 - properties.z_min();

	return	(node.layer() == 1 || node.layer() == 2) &&															// it must be in the lowest branch layer
			MIN_MEAN_GREY_VALUE <= properties.mean_grey_value() && properties.mean_grey_value() <= MAX_MEAN_GREY_VALUE &&	// it must have a reasonable grey value
			properties.x_max() >= volumeSize[0]*0.7 &&										// it must stretch sufficiently far to the right
			properties.y_max() >= spineProperties.y_min() &&										// it should extend as far back as the spine in the y direction			
			properties.voxel_count() >= MIN_VOXELS_PER_SLICE * sliceCount &&					// it must be of a reasonable size
			properties.voxel_count() <= MAX_VOXELS_PER_SLICE * sliceCount;
}

}
