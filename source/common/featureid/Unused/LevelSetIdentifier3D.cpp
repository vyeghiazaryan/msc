/***
 * millipede: LevelSetIdentifier3D.cpp
 * Added by Varduhi Yeghiazaryan, 2013.
 ***/

#include "LevelSetIdentifier3D.h"

#include <boost/bind.hpp>

#include <itkSize.h>

#include <common/dicom/util/WindowSettings.h>
#include <common/dicom/volumes/DICOMVolume.h>
#include <common/interfacemotion/FastMarching.h>

namespace mp {

//#################### CONSTANTS ####################
//const double LevelSetIdentifier3D::LEVELS[] = {1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0};
//const double LevelSetIdentifier3D::LEVELS[] = {2.4e2, 2.1e2, 1.8e2, 1.5e2, 1.2e2, .9e2, .6e2, .3e2, 0.0};
//const double LevelSetIdentifier3D::LEVELS[] = {1.6e1, 1.4e1, 1.2e1, 1.0e1, 0.8e1, 0.6e1, 0.4e1, 0.2e1, 0};
const double LevelSetIdentifier3D::LEVELS[] = {4.8e3, 4.2e3, 3.6e3, 3.0e3, 2.4e3, 1.8e3, 1.2e3, 0.6e3, 0.0};
//const double LevelSetIdentifier3D::LEVELS[] = {8e0, 7e0, 6e0, 5e0, 4e0, 3e0, 2e0, 1e0, 0.0};

//#3################### CONSTRUCTORS ####################
LevelSetIdentifier3D::LevelSetIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF)
:	FeatureIdentifier(dicomVolume, volumeIPF)
{}

//#################### PUBLIC METHODS ####################
int LevelSetIdentifier3D::length() const
{
	return NUMBER + 1;
}

//#################### PRIVATE METHODS ####################
void LevelSetIdentifier3D::execute_impl()
{
	set_status("Identifying level sets...");

	VolumeIPFMultiFeatureSelection_Ptr multiFeatureSelection = get_multi_feature_selection();

	// Step 1: Calculate the combined properties of all the nodes marked as part of the spine.
	BranchProperties spineProperties = multiFeatureSelection->properties_of(AbdominalFeature::VERTEBRA);

	// Step 2: Filter for initial level set.
	std::list<PFNodeID> nodes = filter_branch_nodes(boost::bind(&LevelSetIdentifier3D::is_initial, this, _1, _2, spineProperties));
	std::set<int> leaves;
	for(std::list<PFNodeID>::const_iterator it=nodes.begin(), iend=nodes.end(); it!=iend; ++it)
	{
		std::deque<int> leafIndices = volume_ipf()->receptive_region_of(*it);
		for(std::deque<int>::const_iterator jt=leafIndices.begin(), jend=leafIndices.end(); jt!=jend; ++jt)
		{
			leaves.insert(*jt);
		}	
	}

	// Step 3: Use fast marching method.
	itk::Size<3> volumeSize = dicom_volume()->size();
	std::list<itk::Index<3> > gridIndices;
	for(std::set<int>::const_iterator it=leaves.begin(), iend=leaves.end(); it!=iend; ++it)
	{
		gridIndices.push_back(ITKImageUtil::make_index(	GridUtil::x_of(*it, volumeSize[0]),
								GridUtil::y_of(*it, volumeSize[0], volumeSize[1]),
								GridUtil::z_of(*it, volumeSize[0] * volumeSize[1])	));	
	}
	FastMarching<3, unsigned char> fm(dicom_volume()->windowed_image(WindowSettings(40, 400)), gridIndices);

	// Step 4: Mark the results as level sets.
	for(unsigned int n = 0; n < NUMBER; ++n)
	{
		increment_progress();
		PartitionForestSelection_Ptr region(new PartitionForestSelectionT(volume_ipf()));
		gridIndices = fm.get_shape_at_time(LEVELS[n]);
		for(std::list<itk::Index<3> >::const_iterator it=gridIndices.begin(), iend=gridIndices.end(); it!=iend; ++it)
		{
			int leafIndex = ((*it)[2]*volumeSize[1]+(*it)[1])*volumeSize[0]+(*it)[0];
			region->select_node(PFNodeID(0, leafIndex));
		}
		multiFeatureSelection->identify_selection(region, AbdominalFeature::Enum(AbdominalFeature::LEVEL1+n));
	}
}

bool LevelSetIdentifier3D::is_initial(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties) const
{//Is the same as "is_kidney" from KidneyIdentifier3D!!!
	int sliceCount = properties.z_max() + 1 - properties.z_min();
	int minVoxels = MIN_VOXELS_PER_SLICE * sliceCount;
	int maxVoxels = MAX_VOXELS_PER_SLICE * sliceCount;
	double aspectRatioXY = properties.aspect_ratio_xy();

	return	node.layer() >= 1 &&										// it should be in layer 2 or above (lower nodes aren't kidneys)
			MIN_MEAN_GREY_VALUE <= properties.mean_grey_value() && properties.mean_grey_value() <= MAX_MEAN_GREY_VALUE &&		// it should have a reasonably (but not excessively) high grey value
			(properties.x_max() < spineProperties.centroid().x ||								// it should not cross the spine centroid in the x direction
			 properties.x_min() > spineProperties.centroid().x) &&
			properties.y_max() >= spineProperties.y_min() &&									// it should extend as far back as the spine in the y direction
			MIN_ASPECT_RATIO_XY <= aspectRatioXY && aspectRatioXY <= MAX_ASPECT_RATIO_XY &&										// it should have a reasonable x-y aspect ratio
			minVoxels <= properties.voxel_count() && properties.voxel_count() <= maxVoxels;		// it should be a reasonable size
}

}
