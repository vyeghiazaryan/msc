/***
 * millipede: DICOMCanvas.h
 * Copyright Stuart Golodetz, 2009. All rights reserved.
 ***/

#ifndef H_MILLIPEDE_DICOMCANVAS
#define H_MILLIPEDE_DICOMCANVAS

#include "BaseCanvas.h"

namespace mp {

class DICOMCanvas : public BaseCanvas
{
	//#################### CONSTRUCTORS ####################
public:
	DICOMCanvas(wxWindow *parent, wxGLContext *context, int *attribList, wxWindowID id = -1, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxFULL_REPAINT_ON_RESIZE|wxWANTS_CHARS);

	//#################### PUBLIC METHODS ####################
public:
	void render(wxPaintDC& dc) const;

	//#################### PRIVATE METHODS ####################
private:
	void finish_drawing(wxMouseEvent& e);
	void render_overlays(double left, double top, double right, double bottom) const;
	Greyscale8SliceTextureSet_CPtr texture_set_to_display() const;

	//#################### EVENT HANDLERS ####################
public:
	//~~~~~~~~~~~~~~~~~~~~ MOUSE ~~~~~~~~~~~~~~~~~~~~
	void OnLeaveWindow(wxMouseEvent& e);
	void OnLeftDown(wxMouseEvent& e);
	void OnLeftUp(wxMouseEvent& e);
	void OnMouseMotion(wxMouseEvent& e);
	void OnRightUp(wxMouseEvent& e);

	//#################### EVENT TABLE ####################
	DECLARE_EVENT_TABLE()
};

}

#endif
