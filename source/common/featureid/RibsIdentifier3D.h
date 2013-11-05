/***
 * millipede: RibsIdentifier3D.h
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_RIBSIDENTIFIER3D
#define H_MILLIPEDE_RIBSIDENTIFIER3D

#include <common/jobs/SimpleJob.h>
#include "FeatureIdentifier.h"

namespace mp {

class RibsIdentifier3D : public SimpleJob, public FeatureIdentifier
{
	//#################### CONSTANTS ####################
private:
	static const int MIN_MEAN_GREY_VALUE = 200;
	static const int MIN_VOXELS_PER_SLICE = 90;
	static const int MAX_VOXELS_PER_SLICE = 500;

	//#################### CONSTRUCTORS ####################
public:
	RibsIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF);

	//#################### PUBLIC METHODS ####################
public:
	int length() const;

	//#################### PRIVATE METHODS ####################
private:
	void execute_impl();
	bool grow_condition(const PFNodeID& adj, const BranchProperties& adjProperties, const BranchProperties& curProperties, const BranchProperties& seedProperties, const BranchProperties& overallProperties) const;
	bool is_seed(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties) const;
	PartitionForestSelection_Ptr postprocess_regions(const PartitionForestSelection_Ptr& preliminaryRegions) const;
};

}

#endif
