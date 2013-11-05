/***
 * millipede: AortaIdentifier3D.cpp
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#include "AortaIdentifier3D.h"

#include <cassert>

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>

#include <common/interfacemotion/FastMarching.h>

namespace mp {

//#################### CONSTRUCTORS ####################
AortaIdentifier3D::AortaIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF)
:	FeatureIdentifier(dicomVolume, volumeIPF)
{}

//#################### PUBLIC METHODS ####################
int AortaIdentifier3D::length() const
{
	return 3;
}

//#################### PRIVATE METHODS ####################
void AortaIdentifier3D::execute_impl()
{
	set_status("Identifying aorta...");

	VolumeIPFMultiFeatureSelection_Ptr multiFeatureSelection = get_multi_feature_selection();

	// Step 1: Calculate the combined properties of all the nodes marked as part of the spine and spinal cord.
	BranchProperties spineProperties = multiFeatureSelection->properties_of(AbdominalFeature::VERTEBRA);
	BranchProperties spinalCordProperties = multiFeatureSelection->properties_of(AbdominalFeature::SPINAL_CORD);

	// Step 2: Filter for aorta-like regions.
	std::list<PFNodeID> seeds = filter_branch_nodes(boost::bind(&AortaIdentifier3D::is_seed, this, _1, _2, spineProperties, spinalCordProperties));
	if(seeds.empty()) return;

	PartitionForestSelection_Ptr regions(new PartitionForestSelectionT(volume_ipf()));
	for(std::list<PFNodeID>::const_iterator it=seeds.begin(), iend=seeds.end(); it!=iend; ++it)
	{
		regions->select_node(*it);
	}

	// Step 3: Find the connected components of the results and keep only the one with the greatest yMax.
	int mergeLayer;
	std::set<int> indices;
	boost::tie(mergeLayer, indices) = extract_merge_layer_indices(regions, volume_ipf()->highest_layer());
	std::vector<std::set<int> > connectedComponents = volume_ipf()->find_connected_components(indices, mergeLayer);

	int bestComponentIndex = -1;
	int bestYMax = INT_MIN;
	for(int i=0, size=static_cast<int>(connectedComponents.size()); i<size; ++i)
	{
		BranchProperties componentProperties = calculate_component_properties(mergeLayer, connectedComponents[i]);
		if(componentProperties.y_max() > bestYMax)
		{
			bestComponentIndex = i;
			bestYMax = componentProperties.y_max();
		}
	}

	assert(bestComponentIndex != -1);
	const std::set<int>& bestComponent = connectedComponents[bestComponentIndex];

	// Step 4: Remove any particularly dark regions.
	std::list<PFNodeID> nodes;
	for(std::set<int>::const_iterator it=bestComponent.begin(), iend=bestComponent.end(); it!=iend; ++it)
	{
		PFNodeID node(mergeLayer, *it);
		BranchProperties properties = volume_ipf()->branch_properties(node);
		if(properties.mean_grey_value() >= MIN_MEAN_GREY_VALUE)
		{
			nodes.push_back(node);
		}
	}

	// Step 5: Convert the partition forest candidate nodes to positions in the volume.
	std::list<itk::Index<3> > positions;
	for(std::list<PFNodeID>::const_iterator it=nodes.begin(), iend=nodes.end(); it!=iend; ++it)
	{
		std::deque<int> leaves = volume_ipf()->receptive_region_of(*it);
		for(std::deque<int>::const_iterator jt=leaves.begin(), jend=leaves.end(); jt!=jend; ++jt)
		{
			if(160 < volume_ipf()->leaf_properties(*jt).grey_value() && volume_ipf()->leaf_properties(*jt).grey_value() < 190)
			{
				positions.push_back(volume_ipf()->position_of_leaf(*jt));
			}
		}	
	}

	increment_progress();

	// Step 6: Use the Fast Marching Method.
	FastMarching<3> fm(volume_ipf()->leaf_layer()->gradient_magnitude_image(), positions);
	positions = fm.get_shape_at_time(8.0);

	increment_progress();

	// Step 7: Convert the positions in the volume to partition forest selection.
	std::set<int> leaves;
	for(std::list<itk::Index<3> >::const_iterator it=positions.begin(), iend=positions.end(); it!=iend; ++it)
	{
		leaves.insert(volume_ipf()->leaf_of_position(*it));
	}
	PartitionForestSelection_Ptr region(new PartitionForestSelectionT(volume_ipf(), leaves));

	// Step 8: Mark the results as aorta.
	multiFeatureSelection->identify_selection(region, AbdominalFeature::AORTA);
}

bool AortaIdentifier3D::is_seed(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties,
								const BranchProperties& spinalCordProperties) const
{
	int sliceCount = properties.z_max() + 1 - properties.z_min();
	int minVoxels = MIN_VOXELS_PER_SLICE * sliceCount;
	int maxVoxels = MAX_VOXELS_PER_SLICE * sliceCount;

	return	MIN_MEAN_GREY_VALUE <= properties.mean_grey_value() && properties.mean_grey_value() <= MAX_MEAN_GREY_VALUE &&						// it should have a reasonable grey value
			minVoxels <= properties.voxel_count() && properties.voxel_count() <= maxVoxels &&	// it should have a reasonable size
			spinalCordProperties.x_min() < properties.x_min() && properties.x_min() < spinalCordProperties.x_max() &&						// its left-hand side should be within the spinal cord bounds in x
			properties.y_min() < spineProperties.y_min() &&						// its top should be above the spine
			properties.y_max() > spineProperties.y_min() - 20;					// its bottom should not be much higher than the top of the spine
}

bool AortaIdentifier3D::morphological_condition(const BranchProperties& properties) const
{
	return MIN_MEAN_GREY_VALUE <= properties.mean_grey_value() && properties.mean_grey_value() <= MAX_MEAN_GREY_VALUE;
}

}
