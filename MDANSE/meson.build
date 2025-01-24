project('MDANSE',
        meson_version : '>=1.6')

# List of sources
py_src = {
  'Src/MDANSE/Core': ['__init__.py', 'Error.py', 'Singleton.py', 'SubclassFactory.py',
                      'Platform.py'],
  'Src/MDANSE/MLogging': ['__init__.py'],
  'Src/MDANSE/NeutronInstruments': ['IdealInstrument.py', 'NeutronInstrument.py', '__init__.py'],
  'Src/MDANSE/NeutronInstruments/Resolution': ['Resolution.py', 'IdealResolution.py',
                                               '__init__.py'],
  'Src/MDANSE/NeutronInstruments/Spectrum': ['Spectrum.py', 'FlatSpectrum.py',
                                             '__init__.py'],
  'Src/MDANSE/NeutronInstruments/Method': ['TOFDiffraction.py', 'ScatteringMethod.py',
                                           'NullMethod.py', '__init__.py'],
  'Src/MDANSE/NeutronInstruments/Coverage': [
    'TotalCoverage.py', 'Coverage.py', '__init__.py'
  ],
  'Src/MDANSE/Framework/InputData': ['HDFTrajectoryInputData.py', 'IInputData.py',
                                     'PeriodicTableInputData.py', 'EmptyData.py',
                                     'MVITraceInputData.py', 'InputFileData.py',
                                     '__init__.py', 'MockTrajectoryInputData.py'],
  'Src/MDANSE/Framework/OutputVariables': ['VolumeOutputVariable.py',
                                           'IOutputVariable.py',
                                           'LineOutputVariable.py',
                                           'SurfaceOutputVariable.py',
                                           '__init__.py'],
  'Src/MDANSE/Framework/Projectors': ['PlanarProjector.py', 'IProjector.py',
                                      'NullProjector.py', 'AxialProjector.py',
                                      '__init__.py'],
  'Src/MDANSE/Framework/QVectors': ['LinearQVectors.py',
                                    'MillerIndicesQVectors.py',
                                    'LatticeQVectors.py',
                                    'IQVectors.py',
                                    'CircularLatticeQVectors.py',
                                    'LinearLatticeQVectors.py',
                                    'GridQVectors.py',
                                    'SphericalLatticeQVectors.py',
                                    'DispersionLatticeQVectors.py',
                                    'ApproximateDispersionQVectors.py',
                                    'SphericalQVectors.py',
                                    'CircularQVectors.py',
                                    '__init__.py'],
  'Src/MDANSE/Framework/Formats': ['TextFormat.py', 'IFormat.py', 'SVGFormat.py',
                                   '__init__.py', 'MDTFormat.py', 'MDAFormat.py',
                                   'HDFFormat.py'],
  'Src/MDANSE/Framework/InstrumentResolutions': ['PseudoVoigt.py', 'Triangular.py',
                                                 'Ideal.py', 'Lorentzian.py',
                                                 'Square.py', 'IInstrumentResolution.py',
                                                 '__init__.py', 'Gaussian.py'],
  'Src/MDANSE/Framework/Session': ['CurrentSession.py'],
  'Src/MDANSE/Framework': ['Units.py', 'Configurable.py', 'Status.py', 'UserDefinitionStore.py',
                           '__init__.py'],
  'Src/MDANSE/Framework/Handlers': ['LogfileHandler.py', 'ColorizingStreamHandler.py',
                                    'IHandler.py', '__init__.py'],
  'Src/MDANSE/Framework/Jobs': [
    'GlobalMotionFilteredTrajectory.py', 'AverageStructure.py', 'IJob.py',
    'DipoleAutoCorrelationFunction.py', 'VanHoveFunctionSelf.py', 'UnfoldedTrajectory.py',
    'StaticStructureFactor.py', 'ElasticIncoherentStructureFactor.py', 'JobStatus.py',
    'AreaPerMolecule.py', 'MeanSquareDisplacement.py', 'OrderParameter.py', 'Density.py',
    'DynamicIncoherentStructureFactor.py', 'SolventAccessibleSurface.py',
    'NeutronDynamicTotalStructureFactor.py', 'GaussianDynamicIncoherentStructureFactor.py',
    'CroppedTrajectory.py', 'MolecularTrace.py', 'StructureFactorFromScatteringFunction.py',
    'DensityProfile.py', 'PositionAutoCorrelationFunction.py', 'CoordinationNumber.py',
    'PairDistributionFunction.py', 'Voronoi.py', 'DistanceHistogram.py',
    'PositionPowerSpectrum.py', 'RadiusOfGyration.py', 'Temperature.py',
    'RootMeanSquareFluctuation.py', 'TrajectoryEditor.py', 'Eccentricity.py',
    'MoleculeFinder.py', 'Infrared.py', 'CurrentCorrelationFunction.py',
    'McStasVirtualInstrument.py', 'AngularCorrelation.py', 'DynamicCoherentStructureFactor.py',
    '__init__.py', 'RigidBodyTrajectory.py', 'VelocityAutoCorrelationFunction.py',
    'DensityOfStates.py', 'VanHoveFunctionDistinct.py', 'XRayStaticStructureFactor.py',
    'GeneralAutoCorrelationFunction.py', 'CenterOfMassesTrajectory.py', 'RootMeanSquareDeviation.py'],
  'Src/MDANSE/Framework/AtomSelector': ['molecule_selectors.py', 'all_selector.py',
                                        'group_selectors.py', '__init__.py',
                                        'selector.py', 'atom_selectors.py'],
  'Src/MDANSE/Framework/Converters': ['__init__.py',
    'Forcite.py', 'ASE.py', 'VASP.py', 'MDTraj.py', 'DCD.py',
    'DMol.py', 'Discover.py', 'ImprovedASE.py', 'NAMD.py',
    'CHARMM.py', 'DL_POLY.py', 'MDAnalysis.py', 'Converter.py',
    'CP2K.py', 'LAMMPS.py', 'XPLOR.py', 'DFTB.py', 'Gromacs.py', 'CASTEP.py'],
  'Src/MDANSE/Framework/Configurators': ['__init__.py',
    'AseInputFileConfigurator.py', 'McStasInstrumentConfigurator.py',
    'InputFileConfigurator.py', 'GroupingLevelConfigurator.py',
    'MDTrajTimeStepConfigurator.py', 'XYZFileConfigurator.py',
    'PythonScriptConfigurator.py', 'AtomMappingConfigurator.py',
    'AtomSelectionConfigurator.py', 'MockTrajectoryConfigurator.py',
    'MultiInputFileConfigurator.py', 'MultipleChoicesConfigurator.py',
    'RunningModeConfigurator.py', 'InputDirectoryConfigurator.py',
    'MDAnalysisTimeStepConfigurator.py', 'HDFInputFileConfigurator.py',
    'DistHistCutoffConfigurator.py', 'OutputStructureConfigurator.py',
    'MDFileConfigurator.py', 'TrajectoryVariableConfigurator.py',
    'IntegerConfigurator.py', 'VectorConfigurator.py',
    'FieldFileConfigurator.py', 'OutputFilesConfigurator.py',
    'XTDFileConfigurator.py', 'SingleOutputFileConfigurator.py',
    'McStasOptionsConfigurator.py', 'FramesConfigurator.py',
    'AxisSelectionConfigurator.py', 'BasisSelectionConfigurator.py',
    'InterpolationOrderConfigurator.py', 'QVectorsConfigurator.py',
    'HDFTrajectoryConfigurator.py', 'ProjectionConfigurator.py',
    'MDAnalysisTopologyFileConfigurator.py', 'FileWithAtomDataConfigurator.py',
    'CorrelationFramesConfigurator.py', 'AtomsListConfigurator.py',
    'UnitCellConfigurator.py', 'PythonObjectConfigurator.py',
    'ConfigFileConfigurator.py', 'MoleculeSelectionConfigurator.py',
    'MDTrajTrajectoryFileConfigurator.py', 'AtomTransmutationConfigurator.py',
    'XDATCARFileConfigurator.py', 'DerivativeOrderConfigurator.py',
    'SingleChoiceConfigurator.py', 'FloatConfigurator.py',
    'StringConfigurator.py', 'ASEFileConfigurator.py',
    'McStasParametersConfigurator.py', 'BooleanConfigurator.py',
    'OptionalFloatConfigurator.py', 'PartialChargeConfigurator.py',
    'MDMCTrajectoryConfigurator.py',
    'MDTrajTopologyFileConfigurator.py', 'OutputDirectoryConfigurator.py',
    'MDAnalysisCoordinateFileConfigurator.py', 'OutputTrajectoryConfigurator.py',
    'IConfigurator.py', 'InstrumentResolutionConfigurator.py',
    'RangeConfigurator.py', 'WeightsConfigurator.py'],
  'Src/MDANSE/Framework/AtomMapping': ['__init__.py', 'atom_mapping.py'],
  'Src/MDANSE/IO': ['__init__.py', 'IOUtils.py', 'PDB.py', 'TextFile.py',
                    'FortranFormat.py', 'PDBSpaceGroups.py',
                    'PDBExportFilters.py', 'HDF5.py',
                    'HDF.py', 'PDBReader.py', 'MinimalPDBReader.py'],
  'Src/MDANSE/mdtraj': ['utils.py', '__init__.py'],
  'Src/MDANSE/Extensions': ['__init__.py'],
  'Src/MDANSE/Chemistry': ['Structrures.py', 'Databases.py', 'ChemicalEntity.py',
                           '__init__.py'],
  'Src/MDANSE/MolecularDynamics': ['__init__.py', 'UnitCell.py',
                                   'Configuration.py', 'Analysis.py',
                                   'TrajectoryUtils.py', 'Connectivity.py',
                                   'Trajectory.py', 'MockTrajectory.py',],
  'Src/MDANSE/Scripts': ['mdanse.py', '__init__.py'],
  'Src/MDANSE/Mathematics': ['LinearAlgebra.py', 'Geometry.py', 'Signal.py',
                             '__init__.py', 'Arithmetic.py', 'Transformation.py',
                             'Graph.py'],
  'Src/MDANSE': ['__init__.py'],
  'Src/MDANSE/Trajectory': ['H5MDTrajectory.py', 'MdanseTrajectory.py', '__init__.py']
}

data_files = {
  'Src/MDANSE/Chemistry': ['atoms.json',
                           'molecules.json',
                           'residues_alt_names.json',
                           'nucleotides.json',
                           'residues.json'],
  'Src/MDANSE/Framework': ['units.json']
}

py = import('python').find_installation(pure : false)
dep_py = py.dependency()

np = dependency('numpy', required : true)

foreach folder, sources : py_src
  inc = []
  foreach source : sources
    inc += folder / source
  endforeach
  py.install_sources(inc, subdir : folder.substring(4))  # Strip Src/
endforeach

foreach folder, sources : data_files
  inc = []
  foreach source : sources
    inc += folder / source
  endforeach
  py.install_sources(inc, subdir : folder.substring(4)) # Strip Src/
endforeach

py.install_sources('AUTHORS', 'CHANGELOG', 'CITATION', 'LICENSE', 'README.md',
                   subdir : 'MDANSE')
