//------------------------------------------------------------------------------
//                             OCH TrajectoryData
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// Copyright (c) 2002 - 2022 United States Government as represented by the
// Administrator of The National Aeronautics and Space Administration.
// All Other Rights Reserved.
//
// Author: Jeremy Knittel / NASA
// Created: 2017.02.23
//
/**
 * From original MATLAB prototype:
 *  Author: S. Hughes.  steven.p.hughes@nasa.gov
 */
//------------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include "OCHTrajectoryData.hpp"

//#define DEBUG_READ_OCH
//#define DEBUG_WRITE_OCH

//------------------------------------------------------------------------------
// public methods
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// default constructor
//------------------------------------------------------------------------------
OCHTrajectoryData::OCHTrajectoryData(std::string fileName) :
   TrajectoryData(),
   dataFound(false),
   extraDataFound(false),
   feasibility(0.0),
   optimality(0.0),
   cost(std::numeric_limits<Real>::infinity()),
   maxConViolation(std::numeric_limits<Real>::infinity()),
   exitStatus("N/A"),
   exitCode(0),
   solutionType("N/A")
{
   if (fileName != "")
      ReadFromFile(fileName);
}

//------------------------------------------------------------------------------
// copy constructor
//------------------------------------------------------------------------------
OCHTrajectoryData::OCHTrajectoryData(const OCHTrajectoryData &copy) :
   TrajectoryData(copy),
   dataFound(copy.dataFound),
   extraDataFound(copy.extraDataFound),
   feasibility(copy.feasibility),
   optimality(copy.optimality),
   cost(copy.cost),
   maxConViolation(copy.maxConViolation),
   exitStatus(copy.exitStatus),
   exitCode(copy.exitCode),
   solutionType(copy.solutionType)
{
}

//------------------------------------------------------------------------------
// operator=
//------------------------------------------------------------------------------
OCHTrajectoryData& OCHTrajectoryData::operator=(const OCHTrajectoryData &copy)
{
   
   TrajectoryData::operator=(copy);
   if (&copy == this)
      return *this;      
         
   dataFound = copy.dataFound;
   extraDataFound = copy.extraDataFound;
   feasibility = copy.feasibility;
   optimality = copy.optimality;
   cost = copy.cost;
   maxConViolation = copy.maxConViolation;
   exitStatus = copy.exitStatus;
   exitCode = copy.exitCode;
   solutionType = copy.solutionType;
   return *this;
}

//------------------------------------------------------------------------------
// destructor
//------------------------------------------------------------------------------
OCHTrajectoryData::~OCHTrajectoryData()
{

}

//------------------------------------------------------------------------------
// void SetNumSegments(Integer num)
//------------------------------------------------------------------------------
/**
* This method sets the number of segments in the och data
*
* @param <num> the number of segments
*
* @note We need to override the base definition here since we need the
*       pointers to specifically point to OCHTrajectoryData objects
*
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetNumSegments(Integer num)
{
   bool isLarger = num > numSegments;
   if (isLarger)
   {
      for (Integer ii = 0; ii < (num - numSegments); ii++)
      {
         segments_.push_back(new OCHTrajectorySegment());
         hasSegmentHadDuplicates.push_back(false);
      }
   }
   else
   {
      segments_.resize(num);
      hasSegmentHadDuplicates.resize(num);
   }
   numSegments = num;

}

//------------------------------------------------------------------------------
// void GetCentralBodyAndRefFrame(Integer forSegment,
//                                std::string &itsCentralBody,
//                                std::string &itsRefFrame)
//------------------------------------------------------------------------------
void OCHTrajectoryData::GetCentralBodyAndRefFrame(Integer forSegment,
                                                  std::string &itsCentralBody,
                                                  std::string &itsRefFrame)
{
   if ((forSegment < 0) || (forSegment >= numSegments))
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: segment number ";
      errmsg            += "out of range";
      throw LowThrustException(errmsg);
   }
   itsCentralBody = ((OCHTrajectorySegment*)segments_.at(forSegment))
                    ->GetCentralBody();
   itsRefFrame    = ((OCHTrajectorySegment*)segments_.at(forSegment))
                    ->GetRefFrame();
}

//------------------------------------------------------------------------------
// void SetCentralBodyAndRefFrame(Integer forSegment,
//                                std::string &newCentralBody,
//                                std::string &newRefFrame)
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetCentralBodyAndRefFrame(Integer forSegment,
                                                  std::string newCentralBody,
                                                  std::string newRefFrame)
{
   if ((forSegment < 0) || (forSegment >= numSegments))
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: segment number ";
      errmsg += "out of range";
      throw LowThrustException(errmsg);
   }
   ((OCHTrajectorySegment*)segments_.at(forSegment))
      ->SetCentralBody(newCentralBody);
   ((OCHTrajectorySegment*)segments_.at(forSegment))
      ->SetRefFrame(newRefFrame);
}

//------------------------------------------------------------------------------
// void GetCentralBodyAndRefFrame(Real epoch,
//                                std::string &itsCentralBody,
//                                std::string &itsRefFrame)
//------------------------------------------------------------------------------
void OCHTrajectoryData::GetCentralBodyAndRefFrame(Real epoch,
                                                  std::string &itsCentralBody,
                                                  std::string &itsRefFrame)
{
   Integer relevantSegment = GetRelevantSegment(epoch);
   GetCentralBodyAndRefFrame(relevantSegment, itsCentralBody, itsRefFrame);
}

//------------------------------------------------------------------------------
// void SetTimeSystem(Integer forSegment,
//                    std::string &newTimeSystem)
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetTimeSystem(Integer forSegment,
                                      std::string newTimeSystem)
{
   if ((forSegment < 0) || (forSegment >= numSegments))
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: segment number ";
      errmsg += "out of range";
      throw LowThrustException(errmsg);
   }
   ((OCHTrajectorySegment*)segments_.at(forSegment))
      ->SetTimeSystem(newTimeSystem);
}

//------------------------------------------------------------------------------
// void SetMeshIntervalFractions(Integer forSegment,
//                    Rvector meshIntervalFractions)
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetMeshIntervalFractions(Integer forSegment, 
                                                 Rvector meshIntFracs)
{
   if ((forSegment < 0) || (forSegment >= numSegments))
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: segment number ";
      errmsg += "out of range.";
      throw LowThrustException(errmsg);
   }
   ((OCHTrajectorySegment*)segments_.at(forSegment)) ->
      SetMeshIntervalFractions(meshIntFracs);
}

//------------------------------------------------------------------------------
// std::vector<Rvector> GetAllMeshIntervalFractions()
//------------------------------------------------------------------------------
std::vector<Rvector> OCHTrajectoryData::GetAllMeshIntervalFractions()
{
   std::vector<Rvector> out;
   for (Integer i = 0; i < numSegments; i++)
      out.push_back(((OCHTrajectorySegment*)segments_.at(i))->
         GetMeshIntervalFractions());
   return out;
}

//------------------------------------------------------------------------------
// Rvector GetMeshIntervalFractions(Integer forSegment)
//------------------------------------------------------------------------------
Rvector OCHTrajectoryData::GetMeshIntervalFractions(Integer forSegment)
{
   if ((forSegment < 0) || (forSegment >= numSegments))
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: segment number ";
      errmsg += "out of range.";
      throw LowThrustException(errmsg);
   }
   return ((OCHTrajectorySegment*)segments_.at(forSegment))->
      GetMeshIntervalFractions();
}

//------------------------------------------------------------------------------
// void SetMeshIntervalNumPoints(Integer forSegment,
//                    IntegerArray meshIntervalNumPoints)
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetMeshIntervalNumPoints(Integer forSegment,
                                                 IntegerArray meshIntNumPnts)
{
   if ((forSegment < 0) || (forSegment >= numSegments))
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: segment number ";
      errmsg += "out of range.";
      throw LowThrustException(errmsg);
   }
   ((OCHTrajectorySegment*)segments_.at(forSegment))->
      SetMeshIntervalNumPoints(meshIntNumPnts);
}

//------------------------------------------------------------------------------
// IntegerArray GetMeshIntervalNumPoints(Integer forSegment)
//------------------------------------------------------------------------------
IntegerArray OCHTrajectoryData::GetMeshIntervalNumPoints(Integer forSegment)
{
   if ((forSegment < 0) || (forSegment >= numSegments))
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: segment number ";
      errmsg += "out of range.";
      throw LowThrustException(errmsg);     
   }
   return ((OCHTrajectorySegment*)segments_.at(forSegment))->GetMeshIntervalNumPoints();
}

//------------------------------------------------------------------------------
// std::vector<IntegerArray> GetAllMeshIntervalNumPoints()
//------------------------------------------------------------------------------
std::vector<IntegerArray> OCHTrajectoryData::GetAllMeshIntervalNumPoints()
{
   std::vector<IntegerArray> out;
   for (Integer i = 0; i < numSegments; i++)
      out.push_back(((OCHTrajectorySegment*)segments_.at(i))->
         GetMeshIntervalNumPoints());
   return out;
}

//------------------------------------------------------------------------------
// bool AllMeshDataSet()
//------------------------------------------------------------------------------
bool OCHTrajectoryData::AllMeshDataSet()
{
   bool flag = true;
   for (Integer i = 0; i < numSegments; i++)
      if (!((OCHTrajectorySegment*)segments_.at(i))->MeshDataSet()) {
         flag = false;
         break;
      }
   return flag;
}

//------------------------------------------------------------------------------
// void SetFeasibility(Real feas)
//------------------------------------------------------------------------------
/**
* Set the feasibility of the trajectory data
*
* @param feas The feasibility
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetFeasibility(Real feas)
{
   feasibility = feas;
}

//------------------------------------------------------------------------------
// Real GetFeasibility()
//------------------------------------------------------------------------------
/**
* Get the feasibility of the trajectory data
*
* @return The feasibility
*/
//------------------------------------------------------------------------------
Real OCHTrajectoryData::GetFeasibility()
{
   return feasibility;
}

//------------------------------------------------------------------------------
// void SetOptimality(Real opt)
//------------------------------------------------------------------------------
/**
* Set the optimality of the trajectory data
*
* @param opt The optimality
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetOptimality(Real opt)
{
   optimality = opt;
}

//------------------------------------------------------------------------------
// Real GetOptimality()
//------------------------------------------------------------------------------
/**
* Get the optimality of the trajectory data
*
* @return The optimality
*/
//------------------------------------------------------------------------------
Real OCHTrajectoryData::GetOptimality()
{
   return optimality;
}

//------------------------------------------------------------------------------
// void SetCostValue(Real costVal)
//------------------------------------------------------------------------------
/**
* Set the value of the cost function of the trajectory data
*
* @param costVal The cost value
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetCostValue(Real costVal)
{
   cost = costVal;
}

//------------------------------------------------------------------------------
// Real GetCostValue()
//------------------------------------------------------------------------------
/**
* Set the value of the cost function of the trajectory data
*
* @return The cost value
*/
//------------------------------------------------------------------------------
Real OCHTrajectoryData::GetCostValue()
{
   return cost;
}

//------------------------------------------------------------------------------
// void SetMaxConstraintViolation(Real maxConViol)
//------------------------------------------------------------------------------
/**
* Set the max constraint violation of the trajectory
*
* @param maxConViol The max constraint violation
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetMaxConstraintViolation(Real maxConViol)
{
   maxConViolation = maxConViol;
}

//------------------------------------------------------------------------------
// Real GetMaxConstraintViolation()
//------------------------------------------------------------------------------
/**
* Get the max constraint violation of the trajectory
*
* @return The max constraint violation
*/
//------------------------------------------------------------------------------
Real OCHTrajectoryData::GetMaxConstraintViolation()
{
   return maxConViolation;
}

//------------------------------------------------------------------------------
// void SetOptimizerExitStatus(std::string optExitStatus)
//------------------------------------------------------------------------------
/**
* Set the optimizer exit status that was achieved for the trajectory
*
* @param optExitStatus The exit status
* @param optExitCode The exit code returned from the optimizer
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetOptimizerExitStatus(std::string optExitStatus, 
                                               Integer optExitCode)
{
   exitStatus = optExitStatus;
   exitCode = optExitCode;
}

//------------------------------------------------------------------------------
// std::string GetOptimizerExitStatus()
//------------------------------------------------------------------------------
/**
* Get the optimizer exit status that was achieved for the trajectory
*
* @return The exit status
*/
//------------------------------------------------------------------------------
std::string OCHTrajectoryData::GetOptimizerExitStatus()
{
   return exitStatus;
}

//------------------------------------------------------------------------------
// void SetMaxMeshRelativeErrors(Real maxRelativeMeshErrors)
//------------------------------------------------------------------------------
/**
* Set the max relative mesh error from the trajectory.
*
* @param maxRelMeshErrors The max realative error
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetMaxMeshRelativeError(Real maxRelativeMeshError)
{
   maxRelMeshError = maxRelativeMeshError;
}

//------------------------------------------------------------------------------
// Real GetMaxMeshRelativeErrors()
//------------------------------------------------------------------------------
/**
* Get the max relative mesh error from the trajectory.
*
* @return The max error of each segment
*/
//------------------------------------------------------------------------------
Real OCHTrajectoryData::GetMaxMeshRelativeError()
{
   return maxRelMeshError;
}

//------------------------------------------------------------------------------
// void SetSolutionType(std::string solType)
//------------------------------------------------------------------------------
/**
* Set the solution type that is being stored.
*
* @param <solType> The solution type
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::SetSolutionType(std::string solType)
{
   solutionType = solType;
}

//------------------------------------------------------------------------------
// std::string SetSolutionType()
//------------------------------------------------------------------------------
/**
* Get the solution type that is being stored.
*
* @return The solution type
*/
//------------------------------------------------------------------------------
std::string OCHTrajectoryData::GetSolutionType()
{
   return solutionType;
}

//------------------------------------------------------------------------------
// Real ProcessTimeString(std::string input, std::string timeSystem)
//------------------------------------------------------------------------------
/**
* This method processes a time string based on the time system specified
*
* @param <input>      the time value input
* @param <timeSystem> the time system format of the input
* @return the time value converted to A1ModJulian
*
*/
//------------------------------------------------------------------------------
Real OCHTrajectoryData::ProcessTimeString(std::string input,
                                          std::string timeSystem)
{
   // TODO: process time string and output Real time value.
   // For now just convert it to a string

   Real value;
   bool success = GmatStringUtil::ToReal(input,&value);
   if (success)
   {
      std::string system, format;
      theTimeConverter->
         GetTimeSystemAndFormat(timeSystem, system,
            format);
      Integer origTimeTypeID =
         theTimeConverter->GetTimeTypeID(system);
      std::string a1ModJulainName = "A1";
      Integer a1ModJulianTypeID =
         theTimeConverter->GetTimeTypeID(a1ModJulainName);
      value = theTimeConverter->Convert(value, origTimeTypeID,
         a1ModJulianTypeID);
      return value;
   }
   else
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: error converting ";
      errmsg            += "time string to Real";
      throw LowThrustException(errmsg);
   }
}

//------------------------------------------------------------------------------
// void WriteToFile(std::string fileName)
//------------------------------------------------------------------------------
/**
* This method writes all OCH data to a text file
*
* @param <fileName> the file to write to
*
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::WriteToFile(std::string fileName)
{
   #ifdef DEBUG_WRITE_OCH
      MessageInterface::ShowMessage(" Entering WriteToFile with %s\n",
                                    fileName.c_str());
      MessageInterface::ShowMessage("  number of segments: %d\n",
                                    (Integer) segments_.size());
   #endif
   std::ofstream fOut;
   fOut.open(fileName);
 
   // Get a timestamp to put at the top of the file
   time_t rawtime;
   struct tm * timeinfo;
   time (&rawtime);
   timeinfo = localtime(&rawtime);
   char buffer[80];
   strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);

   // Open the file
   if (fOut.is_open())
   {
      fOut << "Optimal Control History file written by CSALT, "
           << buffer << std::endl;
      fOut << std::endl;

      // Write out the exit conditions of the optimization
      fOut << "EXIT_CONDITIONS_START" << std::endl;
      fOut << "\tOPTIMIZER_FEASIBILITY_METRIC\t= " << std::scientific << std::setprecision(17) << feasibility << std::endl;
      fOut << "\tOPTIMIZER_OPTIMALITY_METRIC\t\t= " << optimality << std::endl;
      fOut << "\tCOST\t\t\t\t\t\t\t= " << cost << std::endl;
      fOut << "\tMAX_CONSTRAINT_VIOLATION\t\t= " << maxConViolation << std::endl;
      fOut << "\tOPTIMIZER_EXIT_STATUS\t\t\t= Exited with code " << exitCode << ", " << exitStatus << std::endl;
      fOut << "\tMAX_RELATIVE_MESH_ERROR\t\t\t= " << maxRelMeshError << std::endl;
      fOut << "\tSOLUTION_TYPE\t\t\t\t\t= " << solutionType << std::endl;
      fOut << "EXIT_CONDITIONS_END" << std::endl;
      fOut << std::endl;
      fOut.unsetf(std::ios::scientific);

      // Loop through all segments (phases)
      for (UnsignedInt s = 0; s < segments_.size(); s++)
      {
         OCHTrajectorySegment *och = (OCHTrajectorySegment*) segments_.at(s);
         #ifdef DEBUG_WRITE_OCH
            MessageInterface::ShowMessage(
                      " Grabbed OCHTrajectorySegment at %d <%p>\n", s, och);
         #endif
        
         // Print the header items
         fOut << "META_START" << std::endl;
         if (och->GetCentralBody().compare(""))  // non-blank central body
            fOut << "\tCENTRAL_BODY\t= " << och->GetCentralBody() << std::endl;
         if (och->GetObjectId().compare(""))     // non-blank ID
            fOut << "\tOBJECT_ID\t\t= " << och->GetObjectId() << std::endl;
         if (och->GetObjectName().compare(""))   // non-blank name
            fOut << "\tOBJECT_NAME\t\t= " << och->GetObjectName() << std::endl;
         if (och->GetRefFrame().compare(""))     // non-blank reference frame
            fOut << "\tREF_FRAME\t\t= " << och->GetRefFrame() << std::endl;
         if (och->GetTimeSystem().compare(""))   // non-blank time system
            fOut << "\tTIME_SYSTEM\t\t= " << och->GetTimeSystem() << std::endl;
         fOut << "\tNUM_STATES\t\t= " << och->GetNumStates() << std::endl;
         fOut << "\tNUM_CONTROLS\t= " << och->GetNumControls() << std::endl;
         fOut << "\tNUM_INTEGRALS\t= " << och->GetNumIntegrals() << std::endl;
         fOut << "META_STOP" << std::endl;
         fOut << std::endl;

         // Print mesh start header
         fOut << "MESH_START" << std::endl;

         // Get mesh interval fractions and number of points per subinterval
         Rvector meshIntervalFractions = och->GetMeshIntervalFractions();
         IntegerArray meshIntervalNumPoints = och->GetMeshIntervalNumPoints();
         Integer num_meshIntFracs = meshIntervalFractions.GetSize();
         Integer num_meshIntNumPnts = meshIntervalNumPoints.size();

         // Loop through mesh interval fractions and print
         for (Integer idx = 0; idx < num_meshIntFracs; idx++)
            if (idx == 0)  // first iteration, begin with tab
               fOut << "\t" << std::setprecision(17) << meshIntervalFractions(idx) << " ";
            else if (idx == num_meshIntFracs - 1)  // last iteration, do not include final space
               fOut << std::setprecision(17) << meshIntervalFractions(idx);
            else 
               fOut << std::setprecision(17) << meshIntervalFractions(idx) << " ";
         fOut << std::endl;

         // Loop through number of points per interval
         for (Integer idx = 0; idx < num_meshIntNumPnts; idx++)
            if (idx == 0)
               fOut << "\t" << meshIntervalNumPoints.at(idx) << " ";
            else if (idx == num_meshIntNumPnts - 1)
               fOut << meshIntervalNumPoints.at(idx);
            else
               fOut << meshIntervalNumPoints.at(idx) << " ";
         fOut << std::endl;

         // Get initial and final time for printing
         Integer numData = och->GetNumDataPoints();
         Real initialTime = och->GetTime(0);
         Real finalTime = och->GetTime(numData - 1);
         fOut << "\t" << std::setprecision(17) << initialTime << " ";
         fOut << std::setprecision(17) << finalTime << std::endl;

         // Finish meta block
         fOut << "MESH_STOP" << std::endl;
         fOut << std::endl;

         // Loop through all data columns and print the state, then control,
         // then integral data
         fOut << "DATA_START" <<std::endl;
         //Integer numData      = och->GetNumDataPoints();
         Integer numStates    = och->GetNumStates();
         Integer numControls  = och->GetNumControls();
         Integer numIntegrals = och->GetNumIntegrals();
         #ifdef DEBUG_WRITE_OCH
            MessageInterface::ShowMessage(
                              " About to loop over numDataPoints\n");
         #endif

         for (Integer idx = 0; idx < numData; idx++)
         {
            Real theTime = och->GetTime(idx);
            fOut << "\t" << std::setw(26) << std::setprecision(17) << std::left
                         << theTime;

            for (Integer stateIdx = 0; stateIdx < numStates; stateIdx++)
            {
               if (och->GetStateSize(idx))   // non-zero state size
                  fOut << " " << std::setw(26) << std::setprecision(17)
                       << std::left << och->GetState(idx,stateIdx);
               else
               {
                  Rvector neededTime;
                  neededTime.SetSize(1);
                  neededTime(0) = theTime;
                  SetAllowExtrapolation(true);
                  SetAllowInterSegmentExtrapolation(true);
                  Rmatrix interpState = GetState(neededTime);
                  fOut << " " << std::setw(26) << std::setprecision(17)
                     << std::left << interpState(0, stateIdx);
               }
            }

            for (Integer controlIdx = 0; controlIdx < numControls; controlIdx++)
            {
               if (och->GetControlSize(idx))
                  fOut << " " << std::setw(26) << std::setprecision(17)
                       << std::left << och->GetControl(idx,controlIdx);
               else
               {
                  Rvector neededTime;
                  neededTime.SetSize(1);
                  neededTime(0) = theTime;
                  SetAllowExtrapolation(true);
                  SetAllowInterSegmentExtrapolation(true);
                  Rmatrix interpControl = GetControl(neededTime);
                  fOut << " " << std::setw(26) << std::setprecision(17)
                       << std::left << interpControl(0,controlIdx);
               }
            }
            
            for (Integer integralIdx = 0; integralIdx < numIntegrals;
                 integralIdx++)
            {
               if (och->GetIntegralSize(idx))
                  fOut << " " << std::setw(26) << std::setprecision(17)
                       << std::left << och->GetIntegral(idx,integralIdx);
               else
               {
                  Rvector neededTime;
                  neededTime.SetSize(1);
                  neededTime(0) = theTime;
                  SetAllowExtrapolation(true);
                  SetAllowInterSegmentExtrapolation(true);
                  Rmatrix interpIntegral = GetIntegral(neededTime);
                  fOut << " " << std::setw(26) << std::setprecision(17)
                       << std::left << interpIntegral(0,integralIdx);
               }
            }
             
            fOut << " " << std::endl;
         } // loop over numData
         fOut << "DATA_STOP" <<std::endl;
         fOut << std::endl;
      } // loop through segments
      fOut.close();
   } // if file is open
   #ifdef DEBUG_WRITE_OCH
   else
   {
      MessageInterface::ShowMessage(" File %s did NOT open!\n",
                                    fileName.c_str());
   }
   #endif
   #ifdef DEBUG_WRITE_OCH
      MessageInterface::ShowMessage(" EXITING WriteToFile with %s\n",
                                    fileName.c_str());
   #endif
}

//------------------------------------------------------------------------------
// void ReadFromFile(std::string fileName)
//------------------------------------------------------------------------------
/**
* This method populates OCH data from an OCH file
*
* @param <fileName> the file to read
*
*/
//------------------------------------------------------------------------------
void OCHTrajectoryData::ReadFromFile(std::string fileName)
{
   #ifdef DEBUG_READ_OCH
      MessageInterface::ShowMessage(" Entering ReadFromFile with %s\n",
                                    fileName.c_str());
   #endif

   // Create and initialize the local variables
   std::ifstream fIn;
   std::string currLine;
   Integer currSegment = -1;
   Integer lineNumber = 0;
   Integer dataLine, meshLine;
   Integer fileLocation = 0; // 1 = in segment header, 2 = in data, 3 = in mesh, 0 = in nothing
   bool    metaFound = false;
   Integer nmStates    = -1;
   Integer nmControls  = -1;
   Integer nmIntegrals = -1;
   TrajectoryDataStructure localData;

   // Clear the vectors
   for (Integer ii = 0; ii < segments_.size(); ii++)
      delete segments_.at(ii);
   segments_.clear();
   hasSegmentHadDuplicates.clear();
   
   numSegments = 0;

   // Open the file
   fIn.open(fileName);
   if (fIn.is_open())
   {
      #ifdef DEBUG_READ_OCH
         MessageInterface::ShowMessage("  ReadFromFile File is open: %s\n",
                                       fileName.c_str());
      #endif
      OCHTrajectorySegment *och = NULL;
      
      // Loop through all lines in the file
      while (std::getline(fIn,currLine))
      {
         lineNumber++;
         #ifdef DEBUG_READ_OCH
            MessageInterface::ShowMessage("  currLine: %s\n", currLine.c_str());
         #endif
        
         if (currLine.find("\"COMMENT") != std::string::npos)
             currLine = currLine.substr(0,currLine.find("\"COMMENT"));            
         
         // Check if a segment header is starting
         if (currLine.find("META_START") != std::string::npos) 
         {
            if ((fileLocation == 1) || // already in META - error
                (fileLocation == 2) || // still in DATA - error
                (fileLocation == 3))   // sill in MESH - error
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Unexpected META_START found in META, ";
               errmsg += "DATA, or MESH block. \n";
               throw LowThrustException(errmsg);
            }
            fileLocation = 1;
            metaFound = true;
            currSegment++;

            // All segments added for this class must be of type
            // OCHTrajectorySegment
            och = new OCHTrajectorySegment();

            // push the new one onto the list at index currSegment
            segments_.push_back(och);
            hasSegmentHadDuplicates.push_back(false);
            numSegments++;
            #ifdef DEBUG_READ_OCH
               MessageInterface::ShowMessage(
                                  "  OCH [%d] <%p> added to the list!!!!!\n",
                                  currSegment, och);
            #endif

         }
         // Check if a segment header is ending
         else if (currLine.find("META_STOP") != std::string::npos)
         {
            fileLocation = 0;
         }
         // Check if a data segment is starting
         else if (currLine.find("DATA_START") != std::string::npos) 
         {
            if (fileLocation == 1)  // still in META - error
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Missing META_STOP. \n";
               throw LowThrustException(errmsg);
            }
            if (fileLocation == 2)   // already in DATA - error
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Unexpected DATA_START found within DATA block. \n";
               throw LowThrustException(errmsg);
            }
            if (fileLocation == 3)
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Unexpected DATA_START found within MESH block. \n";
               throw LowThrustException(errmsg);
            }
            if (!metaFound)
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  No META data found or missing META_START. \n";
               throw LowThrustException(errmsg);
            }
            dataFound = true;
            fileLocation = 2;
            dataLine = -1;
            
            // Make sure we knew how many data columns to expect
            if (nmStates == -1)
            {
               std::string errmsg = "ERROR - Missing NUM_STATES field ";
               errmsg            += "in meta data of ";
               errmsg            += "this file: " + fileName + "\n";
               throw LowThrustException(errmsg);
            }
            if (nmControls == -1)
            {
               std::string errmsg = "ERROR - Missing NUM_CONTROLS field ";
               errmsg            += "in meta data of ";
               errmsg            += "this file: " + fileName + "\n";
               throw LowThrustException(errmsg);
            }
            if (nmIntegrals == -1)
            {
               std::string errmsg = "ERROR - Missing NUM_INTEGRALS field ";
               errmsg            += "in meta data of ";
               errmsg            += "this file: " + fileName + "\n";
               throw LowThrustException(errmsg);
            }
            if ((nmStates == 0) && (nmControls == 0) && (nmIntegrals == 0))
            {
               std::string errmsg = "ERROR: Trajectory data could not be ";
               errmsg            += "read from this file: " + fileName + "\n";
               throw LowThrustException(errmsg);
            }
            localData.states.SetSize(nmStates);
            localData.controls.SetSize(nmControls);
            localData.integrals.SetSize(nmIntegrals);
         }
         // Check if a data segment is ending
         else if (currLine.find("DATA_STOP") != std::string::npos) 
         {
            if (fileLocation != 2)  // error
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Unexpected DATA_STOP found outside DATA block. \n";
               throw LowThrustException(errmsg);
            }
            fileLocation = 0;
//            // reset in case there's another Meta/Data segment/section  TBD
//            nmStates    = -1;
//            nmControls  = -1;
//            nmIntegrals = -1;
         }
         // Check if mesh section is starting
         else if (currLine.find("MESH_START") != std::string::npos)
         {
            if (fileLocation == 1)  // still in META - error
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Missing META_STOP. \n";
               throw LowThrustException(errmsg);
            }
            if (fileLocation == 2)   // already in DATA - error
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Unexpected MESH_START found within DATA block. \n";
               throw LowThrustException(errmsg);
            }
            if (fileLocation == 3)
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Unexpected MESH_START found within MESH block. \n";
               throw LowThrustException(errmsg);
            }
            fileLocation = 3;
            meshLine = -1;
         }
         // Check if a data segment is ending
         else if (currLine.find("MESH_STOP") != std::string::npos) 
         {
            if (fileLocation != 3)  // error
            {
               std::string errmsg = "ERROR reading data from ";
               errmsg += "this file: " + fileName;
               errmsg += ".  Unexpected MESH_STOP found outside MESH block. \n";
               throw LowThrustException(errmsg);
            }
            fileLocation = 0;
         }
         else if (fileLocation == 1) // @todo simplify reading META data
         {
            // Process the header data

            if (currLine.find("CENTRAL_BODY") != std::string::npos)
            {
               std::string input = GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAllBlanks(
                                   currLine.substr(currLine.find("=") + 1,
                                   currLine.length()-currLine.find("=") - 1)),
                                   "\r"),"\n");
               och->SetCentralBody(input);
            }
            else if (currLine.find("OBJECT_ID") != std::string::npos)
            {
               std::string input = GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAllBlanks(
                                   currLine.substr(currLine.find("=") + 1,
                                   currLine.length()-currLine.find("=") - 1)),
                                   "\r"),"\n");
               och->SetObjectId(input);
            }
            else if (currLine.find("OBJECT_NAME") != std::string::npos)
            {
               std::string input = GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAllBlanks(
                                   currLine.substr(currLine.find("=") + 1,
                                   currLine.length()-currLine.find("=") - 1)),
                                   "\r"),"\n");
               och->SetObjectName(input);
            }
            else if (currLine.find("REF_FRAME") != std::string::npos)
            {
               std::string input = GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAllBlanks(
                                   currLine.substr(currLine.find("=") + 1,
                                   currLine.length()-currLine.find("=") - 1)),
                                   "\r"),"\n");
               och->SetRefFrame(input);
            }
            else if (currLine.find("TIME_SYSTEM") != std::string::npos)
            {
               std::string input = GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAll(
                                   GmatStringUtil::RemoveAllBlanks(
                                   currLine.substr(currLine.find("=") + 1,
                                   currLine.length()-currLine.find("=") - 1)),
                                   "\r"),"\n");

               if (theTimeConverter->IsValidTimeSystem(input))
                  och->SetTimeSystem(input);
               else
               {
                  std::string errmsg = "ERROR - OCHTrajectoryData: ";
                  errmsg += "error reading TIME SYSTEM from ";
                  errmsg += "this file: " + fileName + ".  Time system "; 
                  errmsg += "\"" + input + "\" is unknown.\n";
                  throw LowThrustException(errmsg);
               }
            }
            else if (currLine.find("NUM_STATES") != std::string::npos)
            {
               Integer value;
               bool success = GmatStringUtil::ToInteger(
                              GmatStringUtil::RemoveAll(
                              GmatStringUtil::RemoveAll(
                              currLine.substr(currLine.find("=") + 1,
                              currLine.length()-currLine.find("=") - 1),
                              "\r"),"\n"),&value);
               if (success)
               {
                  SetNumStateParams(currSegment,value);
                  nmStates = value;
               }
               else
               {
                  std::string errmsg = "ERROR - OCHTrajectoryData: ";
                  errmsg            += "error reading NUM STATES from ";
                  errmsg            += "this file: " + fileName + "\n";
                  throw LowThrustException(errmsg);
               }
            }
            else if (currLine.find("NUM_CONTROLS") != std::string::npos)
            {
               Integer value;
               bool success = GmatStringUtil::ToInteger(
                              GmatStringUtil::RemoveAll(
                              GmatStringUtil::RemoveAll(
                              currLine.substr(currLine.find("=") + 1,
                              currLine.length()-currLine.find("=") - 1),
                              "\r"),"\n"),&value);
               if (success)
               {
                  SetNumControlParams(currSegment,value);
                  nmControls = value;
               }
               else
               {
                  std::string errmsg = "ERROR - OCHTrajectoryData: ";
                  errmsg            += "error reading NUM CONTROLS from ";
                  errmsg            += "this file: " + fileName + "\n";
                  throw LowThrustException(errmsg);
               }
            }
            else if (currLine.find("NUM_INTEGRALS") != std::string::npos)
            {
               Integer value;
               bool success = GmatStringUtil::ToInteger(
                              GmatStringUtil::RemoveAll(
                              GmatStringUtil::RemoveAll(
                              currLine.substr(currLine.find("=") + 1,
                              currLine.length()-currLine.find("=") - 1),
                              "\r"),"\n"),&value);
               if (success)
               {
                  SetNumIntegralParams(currSegment,value);
                  nmIntegrals = value;
               }
               else
               {
                  std::string errmsg = "ERROR - OCHTrajectoryData: error ";
                  errmsg            += "reading NUM_INTEGRALS from ";
                  errmsg            += "this file: " + fileName + "\n";
                  throw LowThrustException(errmsg);
               }
            }
         }
         else if (fileLocation == 2)
         {
            // Process the data line - ignore BLANK lines
            if (GmatStringUtil::IsBlank(currLine))
               continue;

            dataLine++;
            Integer lineLocation;
            Real value;
            
            // How many data columns do we expect on each line?
            Integer numExpected = nmStates + nmControls + nmIntegrals;

            std::istringstream lineStr;
            lineStr.str(currLine);
            
            //Get the time first - read as a string to use existing
            // ProcessTimeString method
            std::string timeStr;
            lineStr >> timeStr;

            localData.time = ProcessTimeString(timeStr,och->GetTimeSystem());
            
            // larger data array to check for additional data
            Rvector dataVec(numExpected * 2);
            Real dataVal;
            // Try to get all of the data expected
            Integer counter = 0;
            while (lineStr >> dataVal)
            {
               dataVec[counter++] = dataVal;
            }
            if (counter < numExpected)
            {
               std::stringstream errmsg1;
               errmsg1.precision(15);
               errmsg1 << "Error reading this file \"" << fileName
                       << "\": expected " << numExpected << " data columns, "
                       << "but found only " << counter << ".\n";
               throw LowThrustException(errmsg1.str());
            }
            if (counter > numExpected)
            {
               std::stringstream errmsg1;
               errmsg1.precision(15);
               errmsg1 << "Error reading this file \"" << fileName
               << "\": expected only " << numExpected << " data columns, "
               << "but found " << counter << ".\n";
               throw LowThrustException(errmsg1.str());
            }
            Integer idx = 0;
            for (Integer ii = 0; ii < nmStates; ii++)
               localData.states(ii) = dataVec(idx++);
            for (Integer ii = 0; ii < nmControls; ii++)
               localData.controls(ii) = dataVec(idx++);
            for (Integer ii = 0; ii < nmIntegrals; ii++)
               localData.integrals(ii) = dataVec(idx++);
            /// Add the localData to the segment
            try
            {
               och->AddDataPoint(localData);
            }
            catch (LowThrustException &lt)
            {
               // Assume exception thrown is for non-monotonic times
               std::string errmsg = "ERROR initializing guess from ";
               errmsg += "input file \"" + fileName + "\": data points are not ";
               errmsg += "in the correct temporal order.\n";
               throw LowThrustException(errmsg);
            }
         } // process the data line (location == 2)
         else if (fileLocation == 3)
         {
            // Process the data line - ignore BLANK lines
            if (GmatStringUtil::IsBlank(currLine))
               continue;

            meshLine++;

            // Allocate requirements
            std::string valStr;
            if (meshLine == 0)
            {
               // Allocate vector to push mesh fractions to
               std::vector<double> meshFracsVec;

               // Create string stream for current line
               std::istringstream lineStr(currLine);

               // Loop and add mesh fractions to vector 
               while (lineStr >> valStr)
                  meshFracsVec.push_back(stod(valStr));

               // Create Rvector for mesh fracs
               Rvector meshFracs(meshFracsVec.size());
               for (Integer idx = 0; idx < meshFracsVec.size(); idx++)
                  meshFracs(idx) = meshFracsVec.at(idx);

               // Add mesh fractions to OCH data
               SetMeshIntervalFractions(currSegment, meshFracs);
            } 
            else if (meshLine == 1)
            {
               // Allocate vector to push mesh number of points to
               IntegerArray meshNumPoints;

               // Create string stream for current line
               std::istringstream lineStr(currLine);

               // Loop and add mesh num points to vector
               while (lineStr >> valStr)
                  meshNumPoints.push_back(stoi(valStr));

               // Add to mesh interval number of points
               SetMeshIntervalNumPoints(currSegment, meshNumPoints);
            }
            else if (meshLine == 2)
            {
               // Create string stream for current line
               std::istringstream lineStr(currLine);

               // Get initial and final time and set
               lineStr >> valStr;
               och->SetStartTime(stod(valStr));
               lineStr >> valStr;
               och->SetStopTime(stod(valStr));
            }
            else 
            {
               // Incorrect number of lines detected in mesh section
               std::string errmsg = "ERROR initializing guess from ";
               errmsg += "input file \"" + fileName + "\": more than ";
               errmsg += "3 lines detected in MESH data.";
               throw LowThrustException(errmsg);
            }
         }
      } // loop through lines in the file
      if (fileLocation == 1)
      {
         std::string errmsg = "ERROR reading OCH file " + fileName;
         errmsg            += ": META_STOP line ";
         errmsg            += "not found.\n";
         throw LowThrustException(errmsg);
      }
      if (fileLocation == 2)
      {
         std::string errmsg = "ERROR reading OCH file " + fileName;
         errmsg            += ": DATA_STOP line ";
         errmsg            += "not found.\n";
         throw LowThrustException(errmsg);
      }
      if (fileLocation == 3)
      {
         std::string errmsg = "ERROR reading OCH file " + fileName;
         errmsg            += ": MESH_STOP line ";
         errmsg            += "not found.\n";
         throw LowThrustException(errmsg);
      }
      if (!metaFound)
      {
         std::string errmsg = "ERROR reading OCH file " + fileName;
         errmsg            += ": META data ";
         errmsg            += "not found.\n";
         throw LowThrustException(errmsg);
      }
      if (!dataFound)
      {
         std::string errmsg = "ERROR: Trajectory data could not be ";
         errmsg            += "read from this file: " + fileName + "\n";
         throw LowThrustException(errmsg);
      }

   } // is file open
   else 
   {
      std::string errmsg = "ERROR - OCHTrajectoryData: cannot open ";
      errmsg            += "this file: " + fileName + "\n";
      throw LowThrustException(errmsg);
   }
}

//------------------------------------------------------------------------------
// std::vector<TrajectoryDataStructure> Interpolate(Rvector requestedTimes)
//------------------------------------------------------------------------------
/**
* This method performs interpolation at a vector of desired time values
*
* @param <requestedTimes> the time values to poll OCH data at
*
* @return A vector of data structures corresponding to each requested time
*
*/
//------------------------------------------------------------------------------
std::vector<TrajectoryDataStructure> OCHTrajectoryData::Interpolate(
                                     Rvector requestedTimes, DataType type)
{
   #ifdef DEBUG_INTERPOLATE
      MessageInterface::ShowMessage("Entering Interpolate wth type = %d\n",
                                    type);
      MessageInterface::ShowMessage("size of requestedTimes = %d\n",
                                    requestedTimes.GetSize());
   #endif
   
   UpdateInterpolator();

   #ifdef DEBUG_INTERPOLATE
      MessageInterface::ShowMessage("---> Interpolator updated\n");
   #endif

   std::vector<TrajectoryDataStructure> output;
   TrajectoryDataStructure localData;

   localData.states.SetSize(0);
   localData.controls.SetSize(0);
   localData.integrals.SetSize(0);

   // Double check that the time values are consecutive
   Real segTimes[2];
   Real tPrecision = 1.0e-10, relT;

   for (Integer s = 1; s < numSegments; s++)
   {
      segTimes[0] = segments_.at(s-1)->GetTime(segments_.at(s-1)->
                    GetNumDataPoints()-1);
      segTimes[1] = segments_.at(s)->GetTime(0);
      relT = segTimes[0] - segTimes[1];
      if (fabs(segTimes[0]) > 0.1)
         relT /= segTimes[0];

      if (relT > tPrecision && !segmentWarningPrinted)
      {
         std::stringstream msg;
         msg.precision(14);
         msg << "WARNING - TrajectoryData: "
            << "Time inputs between segments are not monotonically "
            << "increasing.  For overlapping times, the first segment "
            << "detected containing the requested time will be used for "
            << "interpolation.  For gaps between segments, interpolation "
            << "will be attempted normally.\n ";
             
         MessageInterface::ShowMessage(msg.str());
         segmentWarningPrinted = true;
      }
   }

   for (Integer idx = 0; idx < requestedTimes.GetSize(); idx++)
   {
      Integer currentSegment = GetRelevantSegment(requestedTimes(idx));
      #ifdef DEBUG_INTERPOLATE
         MessageInterface::ShowMessage(
            "---> relevant segment retrieved for index %d at time %12.10f\n",
            idx, requestedTimes(idx));
      #endif

      // start with empty arrays each time
      localData.states.SetSize(0);
      localData.controls.SetSize(0);
      localData.integrals.SetSize(0);
    
      localData.time = requestedTimes(idx);
      bool success = false;
      if (type == ALL || type == STATE)
      {
         localData.states.SetSize(segments_.at(currentSegment)->GetNumStates());
         for (Integer jdx = 0; jdx < segments_.at(currentSegment)->GetNumStates(); jdx++)
         {
            UpdateInterpPoints(currentSegment,requestedTimes(idx),
                                STATE,jdx);

            if (duplicateTimeFound)
            {
               // In this case, assign the requested data as the known
               // data with the nearest time
               break;
            }

            success = interpolator->Interpolate(requestedTimes(idx),
               &localData.states(jdx));
            if (!success)
            {
               throw LowThrustException("ERROR - TrajectoryData: "
                  "Interpolation of state data failed at time " +
                  GmatStringUtil::ToString(requestedTimes(idx), 14) +
                  " at setgment " +
                  GmatStringUtil::ToString(currentSegment, 1) + ".\n");
            }
         }

         if (duplicateTimeFound)
         {
            if (pointToCopy != -1)
            {
               for (Integer jdx = 0;
                  jdx < segments_.at(currentSegment)->GetNumStates(); jdx++)
               {
                  // Duplicate was found, copy data of closest time to requested
                  // time
                  localData.states(jdx) =
                     segments_.at(currentSegment)->GetState(pointToCopy, jdx);
               }
               // Reset duplicate time data
               duplicateTimeFound = false;
               pointToCopy = -1;
            }
            else
            {
               throw LowThrustException("ERROR - TrajectoryData: "
                  "Interpolation of state data failed at time " +
                  GmatStringUtil::ToString(requestedTimes(idx), 14) +
                  " at setgment " +
                  GmatStringUtil::ToString(currentSegment, 1) + ".\n");
            }
         }
      }

      if (type == ALL || type == CONTROL)
      {
         localData.controls.SetSize(
                            segments_.at(currentSegment)->GetNumControls());
         for (Integer jdx = 0;
              jdx < segments_.at(currentSegment)->GetNumControls(); jdx++)
         {
            UpdateInterpPoints(currentSegment,requestedTimes(idx),
                                CONTROL,jdx);

            if (duplicateTimeFound)
            {
               // In this case, assign the requested data as the known
               // data with the nearest time
               break;
            }

            success = interpolator->Interpolate(requestedTimes(idx),
               &localData.controls(jdx));
            if (!success)
            {
               throw LowThrustException("ERROR - TrajectoryData: "
                  "Interpolation of control data failed at time " +
                  GmatStringUtil::ToString(requestedTimes(idx), 14) +
                  " at setgment " +
                  GmatStringUtil::ToString(currentSegment, 1) + ".\n");
            }
         }

         if (duplicateTimeFound)
         {
            if (pointToCopy != -1)
            {
               for (Integer jdx = 0;
                  jdx < segments_.at(currentSegment)->GetNumControls(); jdx++)
               {
                  // Duplicate was found, copy data of closest time to requested
                  // time
                  localData.controls(jdx) =
                     segments_.at(currentSegment)->GetControl(pointToCopy, jdx);
               }
               // Reset duplicate time data
               duplicateTimeFound = false;
               pointToCopy = -1;
            }
            else
            {
               throw LowThrustException("ERROR - TrajectoryData: "
                  "Interpolation of control data failed at time " +
                  GmatStringUtil::ToString(requestedTimes(idx), 14) +
                  " at setgment " +
                  GmatStringUtil::ToString(currentSegment, 1) + ".\n");
            }
         }
      }

      if (type == ALL || type == INTEGRAL)
      {
         localData.integrals.SetSize(
                   segments_.at(currentSegment)->GetNumIntegrals());
         for (Integer jdx = 0;
              jdx < segments_.at(currentSegment)->GetNumIntegrals(); jdx++)
         {
            UpdateInterpPoints(currentSegment,requestedTimes(idx),
                                INTEGRAL,jdx);

            if (duplicateTimeFound)
            {
               // In this case, assign the requested data as the known
               // data with the nearest time
               break;
            }

            success = interpolator->Interpolate(requestedTimes(idx),
               &localData.integrals(jdx));
            if (!success)
            {
               throw LowThrustException("ERROR - TrajectoryData: "
                  "Interpolation of control data failed at time " +
                  GmatStringUtil::ToString(requestedTimes(idx), 14) +
                  " at setgment " +
                  GmatStringUtil::ToString(currentSegment, 1) + ".\n");
            }
         }

         if (duplicateTimeFound)
         {
            if (pointToCopy != -1)
            {
               for (Integer jdx = 0;
                  jdx < segments_.at(currentSegment)->GetNumIntegrals(); jdx++)
               {
                  // Duplicate was found, copy data of closest time to requested
                  // time
                  localData.integrals(jdx) =
                     segments_.at(currentSegment)->GetIntegral(pointToCopy, jdx);
               }
               // Reset duplicate time data
               duplicateTimeFound = false;
               pointToCopy = -1;
            }
            else
            {
               throw LowThrustException("ERROR - TrajectoryData: "
                  "Interpolation of integral data failed at time " +
                  GmatStringUtil::ToString(requestedTimes(idx), 14) +
                  " at setgment " +
                  GmatStringUtil::ToString(currentSegment, 1) + ".\n");
            }
         }
      }

      output.push_back(localData);
    }

    return output;
}