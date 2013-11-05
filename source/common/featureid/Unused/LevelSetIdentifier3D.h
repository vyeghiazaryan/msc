/***
 * millipede: LevelSetIdentifier3D.h
 * Added by Varduhi Yeghiazaryan, 2013.
 ***/

#ifndef H_MILLIPEDE_LEVELSETIDENTIFIER3D
#define H_MILLIPEDE_LEVELSETIDENTIFIER3D

#include <common/jobs/SimpleJob.h>
#include "FeatureIdentifier.h"

namespace mp {

class LevelSetIdentifier3D : public SimpleJob, public FeatureIdentifier
{
	//#################### CONSTANTS ####################
private:
	static const unsigned short NUMBER = 9;
	static const double LEVELS[NUMBER];

	static const int MIN_MEAN_GREY_VALUE = 150;
	static const int MAX_MEAN_GREY_VALUE = 200;
	static const int MIN_VOXELS_PER_SLICE = 100;
	static const int MAX_VOXELS_PER_SLICE = 3000;
	static const double MIN_ASPECT_RATIO_XY = 0.25;
	static const double MAX_ASPECT_RATIO_XY = 4.0;

	//#################### CONSTRUCTORS ####################
public:
	LevelSetIdentifier3D(const DICOMVolume_CPtr& dicomVolume, const VolumeIPF_Ptr& volumeIPF);

	//#################### PUBLIC METHODS ####################
public:
	int length() const;

	//#################### PRIVATE METHODS ####################
private:
	void execute_impl();
	bool is_initial(const PFNodeID& node, const BranchProperties& properties, const BranchProperties& spineProperties) const;
};

}

#endif
