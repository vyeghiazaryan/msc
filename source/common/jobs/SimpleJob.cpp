/***
 * millipede: SimpleJob.cpp
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 ***/

#include "SimpleJob.h"

namespace mp {

//#################### CONSTRUCTORS ####################
SimpleJob::SimpleJob()
:	m_progress(0)
{}

//#################### PUBLIC METHODS ####################
void SimpleJob::execute()
{
	execute_impl();
	set_progress(length());
}

int SimpleJob::progress() const
{
	boost::mutex::scoped_lock lock(m_mutex);
	return m_progress;
}

std::string SimpleJob::status() const
{
	boost::mutex::scoped_lock lock(m_mutex);
	return m_status;
}

//#################### PROTECTED METHODS ####################
void SimpleJob::increment_progress()
{
	boost::mutex::scoped_lock lock(m_mutex);
	++m_progress;
}

void SimpleJob::set_progress(int progress)
{
	boost::mutex::scoped_lock lock(m_mutex);
	m_progress = progress;
}

}
