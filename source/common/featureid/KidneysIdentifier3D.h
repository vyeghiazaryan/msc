/***
 * millipede: KidneysIdentifier3D.h
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_KIDNEYSIDENTIFIER3D
#define H_MILLIPEDE_KIDNEYSIDENTIFIER3D

#include <common/jobs/SimpleJob.h>
#include "FeatureIdentifier.h"

namespace mp {

class KidneysIdentifier3D : public SimpleJob, public FeatureIdentifier
{
	//#################### CONSTANTS ####################
private:
	static const int MIN_MEAN_GREY_VALUE = 150;
	static const int MAX_MEAN_GREY_VALUE = 218;
	static const int MIN_VOXELS_PER_SLICE = 1800;
	static const int MAX_VOXELS_PER_SLICE = 6000;
	static const double MIN_ASPECT_RATIO_XY = 0.6;
	static const double MAX_ASPECT_RATIO_XY = 1.8;
	static const int MAX_DISTANCE_FROM_SPINE = 150;

	//#################### CONSTRUCTORS ####################
public:
	KidneysIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF);

	//#################### PUBLIC METHODS ####################
public:
	int length() const;

	//#################### PRIVATE METHODS ####################
private:
	void execute_impl();
	bool is_kidney(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties) const;
};

}

#endif
