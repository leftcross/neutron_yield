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
//    For information related to this code contact: M. Settimo
//    e-mail: settimo@subatech.in2p3.fr
// --------------------------------------------------------------
// Comments
//
//                   Underground Advanced -  Modified by M.Settimo 2017
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//               Modified by M.Settimo 2017
//
// DetectorConstruction program
// --------------------------------------------------------------

#include "DAMICDetectorConstruction.hh"
#include "DAMICDetectorMessenger.hh"

#include "DAMICSiSD.hh"
#include "DAMICDetectorConstruction.hh"
#include "DAMICDetectorModules.hh"
#include "G4NistManager.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DAMICDetectorConstruction::DAMICDetectorConstruction()  
{
  // create commands for interactive definition of time cuts:
  detectorMessenger = new DAMICDetectorMessenger(this);

  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theMinEkine         = 250.0*eV; // minimum kinetic energy required in volume

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DAMICDetectorConstruction::~DAMICDetectorConstruction() 
{

  delete detectorMessenger;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DAMICDetectorConstruction::DefineMaterials() 
{

#include "DAMICDetectorMaterial.icc"

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DAMICDetectorConstruction::Construct() {

  DefineMaterials();

  // DefineField();

  // make colours
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.85, .85, .85) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  orange  (.75, .55, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  G4Colour  lgreen  (0.0, .75, 0.0) ;
  G4Colour  green   (0.0, 1.0, 0.0) ;
  G4Colour  brown   (0.7, 0.4, 0.1) ;


  // Universe - room wall - CONCRETE ****************************************

  G4double Offset_DAMICGeoZ = 0 ;  
  G4double wallThick   = 1.*cm;
  G4double worldWidth  = 1600.0*cm + 2.*wallThick; // "x"
  G4double worldLength = 1600.0*cm + 2.*wallThick; // "y"
  G4double worldHeight = 800.0*cm + 2.*wallThick; // "z"


  G4Box* WorldBox = new G4Box("world", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight);
  WorldLV  = new G4LogicalVolume(WorldBox, LabMat, "WorldLV");
  WorldPV = new G4PVPlacement(0, G4ThreeVector(0.,0.,0), "WorldPV", WorldLV, NULL, false,0);

  G4VisAttributes* world_vat= new G4VisAttributes(orange);
  world_vat->SetVisibility(true);
  
  
  // Lab Space - AIR ********************************************************

  G4double labWidth  = worldWidth  - 2.*wallThick; //X
  G4double labLength = worldLength - 2.*wallThick; //Y
  G4double labHeight = worldHeight - 2.*wallThick; //Z

  G4Box* LabBox = new G4Box("Lab_Box", 0.5*labWidth, 0.5*labLength, 0.5*labHeight);
  LabLV = new G4LogicalVolume(LabBox, LabMat, "LabLV");

  G4VisAttributes* lab_vat= new G4VisAttributes(gray);
  lab_vat->SetVisibility(true);
  LabLV->SetVisAttributes(lab_vat);
  

 // Lead box, sorry for naming it as PolyBox, it's actually lead, I didn't change the variable name due to legacy code reason

  G4double OutPolyBoxX = 1580*cm;
  G4double OutPolyBoxY = 1580*cm;
 // G4double OutPolyBoxZ = (16 + 60 + 16)*2.54*cm;
  G4double OutPolyBoxZ = 25*cm;


   G4Box *outerPolyBox = new G4Box("Top face", OutPolyBoxX/2.,OutPolyBoxY/2.,OutPolyBoxZ/2.);
   G4LogicalVolume * extPolyBoxLV = new G4LogicalVolume(outerPolyBox, LeadMat, "extPolyBoxLV");
   G4VisAttributes* poly_vat= new G4VisAttributes(blue);
   poly_vat->SetVisibility(true);
   extPolyBoxLV->SetVisAttributes(poly_vat);


// Position lead box at a gap from origin

  G4double gap=10*cm;   
  G4double outboxcentrez=OutPolyBoxZ/2.+gap;


// Place all the volumes in their corresponding mother volumes to complete the geometry

  LabPV = new G4PVPlacement(0, G4ThreeVector(0.,0.,0), "LabPV", LabLV, WorldPV, false,0);
  G4PVPlacement* extPolyBoxPV = new G4PVPlacement(0, G4ThreeVector(0.,0., outboxcentrez),
						     "extPolyBoxPV", extPolyBoxLV, LabPV, false,true);

  
 
  return WorldPV;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICDetectorConstruction::ConstructSDandField()
{
}
 
void DAMICDetectorConstruction::SetEnergyCut(G4double val)
{

}  



void DAMICDetectorConstruction::SetTimeCut(G4double val)
{
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



