/***
 * millipede: DICOMDIRFile.cpp
 * Copyright Stuart Golodetz, 2009. All rights reserved.
 * Modified by Varduhi Yeghiazaryan, 2013.
 ***/

#include "DICOMDIRFile.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <gdcmReader.h>
#include <gdcmMediaStorage.h>

#include <common/dicom/directories/PatientRecord.h>
#include <common/dicom/directories/SeriesRecord.h>
#include <common/dicom/directories/StudyRecord.h>
#include <common/exceptions/Exception.h>

namespace mp {

//#################### LOADING METHODS ####################
DICOMDirectory_Ptr DICOMDIRFile::load(const std::string& filename)
{
	// Note: The hex values in this function are DICOM tags - they can be looked up in the DICOM standard.

	DICOMDirectory_Ptr ret(new DICOMDirectory);

	typedef std::set<gdcm::DataElement> DataElementSet;
	typedef DataElementSet::const_iterator ConstIterator;

	// Load in the DICOMDIR.
	gdcm::Reader reader;
	reader.SetFileName(filename.c_str());

	if(!reader.Read())
	{
		throw Exception("Failed to load DICOM directory from: " + filename);
	}

	gdcm::File& file = reader.GetFile();
	gdcm::DataSet& ds = file.GetDataSet();
	gdcm::FileMetaInformation& fmi = file.GetHeader();
	
	for(ConstIterator it = ds.GetDES().begin(); it != ds.GetDES().end(); ++it)
		if(it->GetTag() == gdcm::Tag(0x0004, 0x1220))
		{
			const gdcm::DataElement& de = (*it);
			gdcm::SmartPointer<gdcm::SequenceOfItems> sqi = de.GetValueAsSQ();
			unsigned int itemused = 1;
			std::stringstream strm;
			while(itemused <= sqi->GetNumberOfItems())
			{
				strm.str("");
				if (sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0004, 0x1430)))
					sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0004, 0x1430)).GetValue().Print(strm);

				// Add the patients to the directory.
				while((strm.str()=="PATIENT")||((strm.str()=="PATIENT ")))
				{
					strm.str("");
					if(sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0010, 0x0010)))
					sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0010, 0x0010)).GetValue().Print(strm);
					PatientRecord_Ptr patientRecord(new PatientRecord(strm.str()));

					/*ADD TAG TO READ HERE*/
					itemused++;
					strm.str("");
					if(sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0004, 0x1430)))
					sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0004, 0x1430)).GetValue().Print(strm);
					
					// Add the studies to the patient record.
					while((strm.str()=="STUDY")||((strm.str()=="STUDY ")))
					{
						//STUDY DESCRIPTION
						strm.str("");
						if(sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0008, 0x1030)))
						sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0008, 0x1030)).GetValue().Print(strm);
						std::string studyDescription = strm.str();
						
						//STUDY ID
						strm.str("");
						if(sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0020, 0x0010)))
						sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0020, 0x0010)).GetValue().Print(strm);
						std::string studyID = strm.str();
						
						//UID
						strm.str("");
						if (sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0020, 0x000d)))
						sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0020, 0x000d)).GetValue().Print(strm);
						std::string studyInstanceUID = strm.str();

						// Sanitize values (bleurgh - this really shouldn't be necessary).
						studyInstanceUID = studyInstanceUID.substr(0, studyInstanceUID.find_first_of('\0'));

						StudyRecord_Ptr studyRecord(new StudyRecord(studyDescription, studyID, studyInstanceUID));            

						/*ADD TAG TO READ HERE*/
						itemused++;
						strm.str("");
						if(sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0004, 0x1430)))
						sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0004, 0x1430)).GetValue().Print(strm);

						// Add the series to the study record.
						while((strm.str()=="SERIES")||((strm.str()=="SERIES ")))
						{
							//SERIES NUMBER
							strm.str("");
							if (sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0020, 0x0011)))
							sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0020, 0x0011)).GetValue().Print(strm);
							std::string seriesNumber = strm.str();
							boost::trim(seriesNumber);
							SeriesRecord_Ptr seriesRecord(new SeriesRecord(seriesNumber));
							
							/*ADD TAG TO READ HERE*/
							itemused++;
							strm.str("");
							if(sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0004, 0x1430)))
							sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0004, 0x1430)).GetValue().Print(strm);
							
							// Add the image filenames to the series record.
							while ((strm.str()=="IMAGE")||((strm.str()=="IMAGE ")))
							{
								strm.str("");
								if(sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0004, 0x1500)))
								sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0004, 0x1500)).GetValue().Print(strm);
								std::string referencedFileID = strm.str();
								boost::replace_all(referencedFileID, "\\", "/");
								boost::trim(referencedFileID);
								seriesRecord->add_image_filename(referencedFileID);
								
								/*ADD TAG TO READ HERE*/
								if(itemused < sqi->GetNumberOfItems())
									itemused++;
								else break;
								strm.str("");
								if(sqi->GetItem(itemused).FindDataElement(gdcm::Tag (0x0004, 0x1430)))
								sqi->GetItem(itemused).GetDataElement(gdcm::Tag (0x0004, 0x1430)).GetValue().Print(strm);
							}
							studyRecord->add_series_record(seriesRecord);
						}
						patientRecord->add_study_record(studyRecord);
					}
					ret->add_patient_record(patientRecord);
				}			
				++itemused;
			}
		}
/*	// Add the patients to the directory.
	for(gdcm::DicomDirPatient *patient=dicomdir.GetFirstPatient(); patient!=NULL; patient=dicomdir.GetNextPatient())
	{
		std::string patientsName = patient->GetEntryValue(0x0010, 0x0010);
		PatientRecord_Ptr patientRecord(new PatientRecord(patientsName));

		// Add the studies to the patient record.
		for(gdcm::DicomDirStudy *study=patient->GetFirstStudy(); study!=NULL; study=patient->GetNextStudy())
		{
			std::string studyDescription = study->GetEntryValue(0x0008, 0x1030);
			std::string studyID = study->GetEntryValue(0x0020, 0x0010);
			std::string studyInstanceUID = study->GetEntryValue(0x0020, 0x000d);

			// Sanitize values (bleurgh - this really shouldn't be necessary).
			studyInstanceUID = studyInstanceUID.substr(0, studyInstanceUID.find_first_of('\0'));

			StudyRecord_Ptr studyRecord(new StudyRecord(studyDescription, studyID, studyInstanceUID));

			// Add the series to the study record.
			for(gdcm::DicomDirSerie *serie=study->GetFirstSerie(); serie!=NULL; serie=study->GetNextSerie())
			{
				std::string seriesNumber = serie->GetEntryValue(0x0020, 0x0011);
				boost::trim(seriesNumber);
				SeriesRecord_Ptr seriesRecord(new SeriesRecord(seriesNumber));

				// Add the image filenames to the series record.
				for(gdcm::DicomDirImage *image=serie->GetFirstImage(); image!=NULL; image=serie->GetNextImage())
				{
					std::string referencedFileID = image->GetEntryValue(0x0004, 0x1500);
					boost::replace_all(referencedFileID, "\\", "/");
					boost::trim(referencedFileID);
					seriesRecord->add_image_filename(referencedFileID);
				}

				studyRecord->add_series_record(seriesRecord);
			}

			patientRecord->add_study_record(studyRecord);
		}

		ret->add_patient_record(patientRecord);
	}*/

	return ret;
}

}
