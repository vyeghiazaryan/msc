/***
 * millipede: Greyscale8ImageTexture.h
 * Copyright Stuart Golodetz, 2009. All rights reserved.
 ***/

#ifndef H_MILLIPEDE_GREYSCALE8IMAGETEXTURE
#define H_MILLIPEDE_GREYSCALE8IMAGETEXTURE

#include "ITKImageTexture.h"

namespace mp {

//#################### TYPEDEFS ####################
typedef itk::Image<unsigned char,2> Greyscale8Image;

class Greyscale8ImageTexture : public ITKImageTexture<Greyscale8Image>
{
	//#################### FRIENDS ####################
	friend class TextureFactory;

	//#################### CONSTRUCTORS ####################
private:
	Greyscale8ImageTexture(const ImagePointer& image, bool clamp);

	//#################### PUBLIC METHODS ####################
public:
	boost::shared_ptr<ITKImageTexture<Greyscale8Image> > clone() const;

	//#################### PRIVATE METHODS ####################
private:
	void reload_image() const;
	void reload_partial_image(int minX, int minY, int maxX, int maxY) const;
};

//#################### TYPEDEFS ####################
typedef boost::shared_ptr<Greyscale8ImageTexture> Greyscale8ImageTexture_Ptr;

}

#endif
