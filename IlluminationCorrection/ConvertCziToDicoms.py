import os, sys, glob, stat, getopt, platform, subprocess, multiprocessing
'''import Tkinter, Tkconstants, tkFileDialog

class TkFileDialog(Tkinter.Frame):
  def __init__(self, root):

    pform = platform.system()
    self.execExt = ''
    self.execPref = ''
    if pform.lower() == 'windows':
      self.execExt = '.exe'

    Tkinter.Frame.__init__(self, root)

    # options for buttons
    button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}

    # define buttons
    Tkinter.Button(self, text='Shade Correct Set', command=self.shadecorrect).pack(**button_opt)
    Tkinter.Button(self, text='Register Set', command=self.register).pack(**button_opt)
    Tkinter.Button(self, text='Shade Correct+Register Set',
                  command=self.shadecorrectnregister).pack(**button_opt)

    # define options for opening or saving a file
    self.file_opt = options = {}
    options['defaultextension'] = '.czi'
    options['filetypes'] = [('czi files', '.czi'),('all files', '.*')]
    options['parent'] = root
    options['title'] = 'Select files for processing'

    # if you use the multiple file version of the module functions this option is set automatically.
    options['multiple'] = 1

  def askopenfilename( self ):
    # get filename
    filename = tkFileDialog.askopenfilename(**self.file_opt)
    return filename

  def shadecorrect(self):
    files = self.askopenfilename()
    for filename in files:
	self.WriteNrrdTilesFromCziNShadeCorrect(filename)'''

class PassFileNamesFromCLInsteadOfTk():
  def __init__(self, files):
    pform = platform.system()
    self.execExt = ''
    self.execPref = ''
    if pform.lower() == 'windows':
      self.execExt = '.exe'
    else:
      self.execPref = './'
    numCores = multiprocessing.cpu_count()
    self.numCoresToUse = str(int( round(numCores*0.9,0) ));
    skip = 1							#No gui code
    for filename in files:					#
      if skip:							#
        skip = 0						#
      else:							#
        self.WriteNrrdTilesFromCziNShadeCorrect(filename)	#

  def register( self ):
    #files = self.askopenfilename()
    #for filename in files:
    #    print filename 
    print('Nothing implemented yet. Empty for now.')

  def shadecorrectnregister( self ):
    #files = self.askopenfilename()
    #for filename in files:
    #    print filename 
    print('Nothing implemented yet. Empty for now.')

  def WriteNrrdTilesFromCziNShadeCorrect( self, filename ):
    #Get Metadata
    self.RunBioformatsMetaReader(filename)
    #Nrrd files before conversion
    searchStr = os.path.dirname(filename)+'/*.nrrd'
    filesThatExistPre = glob.glob(searchStr)
    #Write Nrrd files
    nrrdConverter = self.execPref+'CziToNrrd'+self.execExt
    args = [ nrrdConverter, filename, self.numCoresToUse ]
    self.RunExec( args, filename )
    #Search and stitch new nrrd files into dicoms
    IlluminationEx = self.execPref+'Illumination'+self.execExt
    dicomConverter = self.execPref+'NrrdToDicom'+self.execExt
    filesThatExistPost = glob.glob(searchStr)
    for nrrdFile in filesThatExistPost:
      noSkip = 1
      if any(nrrdFile in nrrdFileNoConvert for nrrdFileNoConvert in filesThatExistPre):
        noSkip = 0
      if noSkip:
        args = [ IlluminationEx, nrrdFile, self.numCoresToUse ]
#        self.RunExec( args, nrrdFile )
	nrrdIllFile = os.path.splitext(nrrdFile)[0]+'_IlluminationCorrected.nrrd'
        args = [ dicomConverter, nrrdIllFile, os.path.splitext(filename)[0]+'.xml', self.numCoresToUse ]
#        self.RunExec( args, nrrdIllFile )
	print nrrdFile, nrrdIllFile

  def RunBioformatsMetaReader( self, filename ):
    xmlFilename = os.path.splitext(filename)[0]+'.xml'
    errFilename = os.path.splitext(filename)[0]+'.errlog'
    javaExec = 'java'+self.execExt
    args = [ javaExec, '-mx512m', '-cp', 'loci_tools.jar', 'loci.formats.tools.ImageInfo',
           '-nopix', '-omexml-only', filename ]
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    f = open( xmlFilename, 'w' )
    f.write(stdout)
    f.close()
    f0 = open( errFilename, 'a' )
    f0.write( stderr )
    f0.close()

  def RunExec( self, args, filename ):
    errFilename = os.path.splitext(filename)[0]+'.errlog'
    logFilename = os.path.splitext(filename)[0]+'.log'
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    f1 = open( errFilename, 'a' )
    f1.write( stderr )
    f1.close()
    f2 = open( logFilename, 'a' )
    f2.write( stdout )
    f2.close()

if __name__=='__main__':
  PassFileNamesFromCLInsteadOfTk(sys.argv)	#No gui code
'''  root = Tkinter.Tk()
  TkFileDialog(root).pack()
  root.mainloop()'''

