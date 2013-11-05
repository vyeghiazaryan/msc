/***
 * millipede: AortaIdentifier3D.h
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_AORTAIDENTIFIER3D
#define H_MILLIPEDE_AORTAIDENTIFIER3D

#include <common/jobs/SimpleJob.h>
#include "FeatureIdentifier.h"

namespace mp {

class AortaIdentifier3D : public SimpleJob, public FeatureIdentifier
{
	//#################### CONSTANTS ####################
private:
	static const int MIN_MEAN_GREY_VALUE = 150;
	static const int MAX_MEAN_GREY_VALUE = 202;
	static const int MIN_VOXELS_PER_SLICE = 350;
	static const int MAX_VOXELS_PER_SLICE = 900;

	//#################### CONSTRUCTORS ####################
public:
	AortaIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF);

	//#################### PUBLIC METHODS ####################
public:
	int length() const;

	//#################### PRIVATE METHODS ####################
private:
	void execute_impl();
	bool is_seed(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties, const BranchProperties& spinalCordProperties) const;
	bool morphological_condition(const BranchProperties& properties) const;
};

}

#endif
