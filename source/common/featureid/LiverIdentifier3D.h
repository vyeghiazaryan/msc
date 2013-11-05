/***
 * millipede: LiverIdentifier3D.h
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_LIVERIDENTIFIER3D
#define H_MILLIPEDE_LIVERIDENTIFIER3D

#include <common/jobs/SimpleJob.h>
#include "FeatureIdentifier.h"

namespace mp {

class LiverIdentifier3D : public SimpleJob, public FeatureIdentifier
{
	//#################### CONSTANTS ####################
private:
	static const int MIN_MEAN_GREY_VALUE = 145;
	static const int MAX_MEAN_GREY_VALUE = 190;
	static const int MIN_VOXELS_PER_SLICE = 700;
	
	//#################### CONSTRUCTORS ####################
public:
	LiverIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF);

	//#################### PUBLIC METHODS ####################
public:
	int length() const;

	//#################### PRIVATE METHODS ####################
private:
	void execute_impl();
	bool grow_condition(const PFNodeID& adj, const BranchProperties& adjProperties, const BranchProperties& curProperties, const BranchProperties& seedProperties, const BranchProperties& overallProperties) const;
	bool is_candidate(const PFNodeID& node, const BranchProperties& properties) const;
	bool morphological_condition(const BranchProperties& properties) const;
};

}

#endif
