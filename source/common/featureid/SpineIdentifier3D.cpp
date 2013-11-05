/***
 * millipede: SpineIdentifier3D.cpp
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#include "SpineIdentifier3D.h"

#include <boost/bind.hpp>

#include <common/dicom/volumes/DICOMVolume.h>
#include <common/interfacemotion/FastMarching.h>
#include <common/util/ITKImageUtil.h>

namespace mp {

//#################### CONSTRUCTORS ####################
SpineIdentifier3D::SpineIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF)
:	FeatureIdentifier(dicomVolume, volumeIPF)
{}

//#################### PUBLIC METHODS ####################
int SpineIdentifier3D::length() const
{
	return 3;
}

//#################### PRIVATE METHODS ####################
void SpineIdentifier3D::execute_impl()
{
	set_status("Identifying the spine...");

	VolumeIPFMultiFeatureSelection_Ptr multiFeatureSelection = get_multi_feature_selection();

	// Step 1: Filter for spine.
	std::list<PFNodeID> nodes = filter_branch_nodes(boost::bind(&SpineIdentifier3D::is_spine, this, _1, _2));

	// Step 2: Convert the partition forest candidate nodes to positions in the volume.
	std::set<int> leaves;
	for(std::list<PFNodeID>::const_iterator it=nodes.begin(), iend=nodes.end(); it!=iend; ++it)
	{
		std::deque<int> leafIndices = volume_ipf()->receptive_region_of(*it);
		for(std::deque<int>::const_iterator jt=leafIndices.begin(), jend=leafIndices.end(); jt!=jend; ++jt)
		{
			leaves.insert(*jt);
		}	
	}
	std::list<itk::Index<3> > positions;
	for(std::set<int>::const_iterator it=leaves.begin(), iend=leaves.end(); it!=iend; ++it)
	{
		if(210 < volume_ipf()->leaf_properties(*it).grey_value())
		{
			positions.push_back(volume_ipf()->position_of_leaf(*it));
		}
	}

	increment_progress();

	// Step 3: Use the Fast Marching Method.
	FastMarching<3> fm(volume_ipf()->leaf_layer()->gradient_magnitude_image(), positions);
	positions = fm.get_shape_at_time(100.0);

	increment_progress();

	// Step 4: Convert the positions in the volume to partition forest selection.
	leaves.clear();
	for(std::list<itk::Index<3> >::const_iterator it=positions.begin(), iend=positions.end(); it!=iend; ++it)
	{
		leaves.insert(volume_ipf()->leaf_of_position(*it));
	}
	PartitionForestSelection_Ptr region(new PartitionForestSelectionT(volume_ipf(), leaves));

	// Step 5: Mark the results as spine.
	multiFeatureSelection->identify_selection(region, AbdominalFeature::VERTEBRA);
}

bool SpineIdentifier3D::is_spine(const PFNodeID& node, const BranchProperties& properties) const
{
	itk::Index<3> volumeSize = ITKImageUtil::make_index_from_size(dicom_volume()->size());
	int minVoxels = MIN_VOXELS_PER_SLICE * volumeSize[2];
	int maxVoxels = MAX_VOXELS_PER_SLICE * volumeSize[2];
	double aspectRatioXY = properties.aspect_ratio_xy();

	return	properties.x_min() < volumeSize[0]/2 && properties.x_max() > volumeSize[0]/2 &&		// it should straddle x = volumeSize[0] / 2
			properties.centroid().y > volumeSize[1]/2 &&												// its centroid should be below y = volumeSize[1]/2
			properties.z_min() == 0 && properties.z_max() == volumeSize[2]-1 &&					// it should extend through all the slices we're looking at
			MIN_ASPECT_RATIO_XY <= aspectRatioXY && aspectRatioXY <= MAX_ASPECT_RATIO_XY &&										// it should have a reasonable x-y aspect ratio
			properties.mean_grey_value() >= MIN_MEAN_GREY_VALUE &&												// it should have a reasonably high grey value
			properties.voxel_count() >= minVoxels && properties.voxel_count() <= maxVoxels;		// and a reasonable size
}

}
