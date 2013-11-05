/***
 * millipede: SpinalCordIdentifier3D.h
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_SPINALCORDIDENTIFIER3D
#define H_MILLIPEDE_SPINALCORDIDENTIFIER3D

#include <common/jobs/SimpleJob.h>
#include "FeatureIdentifier.h"

namespace mp {

class SpinalCordIdentifier3D : public SimpleJob, public FeatureIdentifier
{
	//#################### CONSTANTS ####################
private:
	static const int MAX_MEAN_GREY_VALUE = 140;
	static const int MIN_VOXELS_PER_SLICE = 300;
	static const int MAX_VOXELS_PER_SLICE = 1000;	

	//#################### CONSTRUCTORS ####################
public:
	SpinalCordIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF);

	//#################### PUBLIC METHODS ####################
public:
	int length() const;

	//#################### PRIVATE METHODS ####################
private:
	void execute_impl();
	bool is_spinal_cord(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties) const;
};

}

#endif
