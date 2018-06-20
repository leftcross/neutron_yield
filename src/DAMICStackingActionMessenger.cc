//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Muon Induced Background
//
//      For information related to this code contact: M. Settimo
//    e-mail: settimo@subatech.in2p3.fr
// --------------------------------------------------------------
// Comments
//
//                   Underground Advanced -  Modified by M.Settimo 2017
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// StackingActionMessenger program
// --------------------------------------------------------------

#include "DAMICStackingActionMessenger.hh"

#include <sstream>

#include "DAMICStackingAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"
#include "globals.hh"


DAMICStackingActionMessenger::DAMICStackingActionMessenger
  (DAMICStackingAction* msa) :DMXAction(msa) {

  KillGammasCmd = new G4UIcmdWithABool("/damic/KillGammasInConcrete",this);
  KillGammasCmd->
    SetGuidance("Kills gammas produced by neutrons in the concrete wall");
  KillGammasCmd->SetGuidance("Default = false");
  KillGammasCmd->SetParameterName("KillGammasFlag", false);
  KillGammasCmd->SetDefaultValue(false);

}

DAMICStackingActionMessenger::~DAMICStackingActionMessenger() {

  delete KillGammasCmd;
}

void DAMICStackingActionMessenger::SetNewValue
   (G4UIcommand* command, G4String newValue) {

  if(command==KillGammasCmd) {
    G4int vl;
    const char* t = newValue;
    std::istringstream is(t);
    is >> vl;
    DMXAction->SetKillGammasFlag(vl!=0);
  }

}

