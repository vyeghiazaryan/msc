/***
 * millipede: KidneysIdentifier3D.cpp
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#include "KidneysIdentifier3D.h"

#include <boost/bind.hpp>

#include <common/dicom/volumes/DICOMVolume.h>
#include <common/interfacemotion/FastMarching.h>

namespace mp {

//#################### CONSTRUCTORS ####################
KidneysIdentifier3D::KidneysIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF)
:	FeatureIdentifier(dicomVolume, volumeIPF)
{}

//#################### PUBLIC METHODS ####################
int KidneysIdentifier3D::length() const
{
	return 3;
}

//#################### PRIVATE METHODS ####################
void KidneysIdentifier3D::execute_impl()
{
	set_status("Identifying kidneys...");

	VolumeIPFMultiFeatureSelection_Ptr multiFeatureSelection = get_multi_feature_selection();

	// Step 1: Calculate the combined properties of all the nodes marked as part of the spine.
	BranchProperties spineProperties = multiFeatureSelection->properties_of(AbdominalFeature::VERTEBRA);

	// Step 2: Filter for kidneys.
	std::list<PFNodeID> nodes = filter_branch_nodes(boost::bind(&KidneysIdentifier3D::is_kidney, this, _1, _2, spineProperties));
	
	// Step 3: Convert the partition forest candidate nodes to positions in the volume.
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
		if(160 < volume_ipf()->leaf_properties(*it).grey_value() && volume_ipf()->leaf_properties(*it).grey_value() < 190)
		{
			positions.push_back(volume_ipf()->position_of_leaf(*it));
		}
	}

	increment_progress();

	// Step 4: Use the Fast Marching Method.
	FastMarching<3> fm(volume_ipf()->leaf_layer()->gradient_magnitude_image(), positions);
	int slices = dicom_volume()->size()[2];
	positions = fm.get_shape_at_first_stop(40 * slices);

	increment_progress();

	// Step 5: Convert the positions in the volume to partition forest selection.
	leaves.clear();
	for(std::list<itk::Index<3> >::const_iterator it=positions.begin(), iend=positions.end(); it!=iend; ++it)
	{
		leaves.insert(volume_ipf()->leaf_of_position(*it));
	}
	PartitionForestSelection_Ptr region(new PartitionForestSelectionT(volume_ipf(), leaves));

	// Step 6: Mark the results as kidney.
	multiFeatureSelection->identify_selection(region, AbdominalFeature::KIDNEY);
}

bool KidneysIdentifier3D::is_kidney(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties) const
{
	int sliceCount = properties.z_max() + 1 - properties.z_min();
	int minVoxels = MIN_VOXELS_PER_SLICE * sliceCount;
	int maxVoxels = MAX_VOXELS_PER_SLICE * sliceCount;
	double aspectRatioXY = properties.aspect_ratio_xy();

	return	node.layer() >= 3 &&										// it should be in layer 3 or above (lower nodes aren't kidneys)
			MIN_MEAN_GREY_VALUE <= properties.mean_grey_value() && properties.mean_grey_value() <= MAX_MEAN_GREY_VALUE &&														// it should have a reasonably (but not excessively) high grey value
			(properties.x_max() < spineProperties.centroid().x ||					// it should not cross the spine centroid in the x direction
			 properties.x_min() > spineProperties.centroid().x) &&
			(properties.x_min() > spineProperties.centroid().x - MAX_DISTANCE_FROM_SPINE &&		// it should be close enough to the spine in the x direction
			 properties.x_max() < spineProperties.centroid().x + MAX_DISTANCE_FROM_SPINE) &&
			properties.y_max() >= spineProperties.y_min() &&					// it should extend as far back as the spine in the y direction
			MIN_ASPECT_RATIO_XY <= aspectRatioXY && aspectRatioXY <= MAX_ASPECT_RATIO_XY &&		// it should have a reasonable x-y aspect ratio
			minVoxels <= properties.voxel_count() && properties.voxel_count() <= maxVoxels;		// it should be a reasonable size
}

}
