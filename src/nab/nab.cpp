#include <stdlib.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <cstring>
#include <vector>
#include <sstream>
#include <cstdio>
#include <iostream>

#ifdef WIN32
#	include <fcntl.h>
#	include <windows.h>

#	define STDERR_FILENO 2
#	define STDIN_FILENO 0
#	define STDOUT_FILENO 1
#else
#	include <unistd.h>
#   include <sys/wait.h>
#	include <errno.h>
#endif

#define CPP "ucpp"
#define CREATED_OBJECT_SUFFIX ".o"

//initialized by code generated and linked by the CMake build system
extern std::vector<std::string> builtinLibraries;
extern std::vector<std::string> builtinLinkDirs;
extern std::vector<std::string> builtinCompilerFlags;

static void compileNab(std::vector<std::string> inputFileBasenames, std::vector<std::string> cppArgs);
static void compileC(std::string c_compiler, std::vector<std::string> inputFileBasenames);
static void linkC(std::string linker, std::vector<std::string> objectFiles, std::string outputFilename);

bool compileOnly = false;
bool verbose = false;
bool nab2cDebugMode = false;
bool noassert = false;
bool nodebug = false;
bool saveIntermediates = false;

std::string amberhome;

// Class to handle building command lines for both execvp and CreateProcess
// the main issue is that CreateProcess requires paths with spaces to be quoted, while execvp doesn't understand that syntax.
class CommandLine
{
	std::string _program;

	std::vector<std::string> _arguments;

	std::string quoteIfNecessary(std::string toQuote)
	{
		if(toQuote.find(" ") != std::string::npos)
		{
			toQuote = '\"' + toQuote + '\"';
		}

		return toQuote;
	}

	public:

	//construct with full path to or name of program to execute
	CommandLine(std::string program):
	_program(program),
	_arguments()
	{

	}

	// adds an argument.  This is NOT a simple string concatenation; the argument should be one element of the target program's argv array.
	// Spaces will be quoted automatically if necessary.
	CommandLine& arg(std::string arg)
	{
		_arguments.push_back(arg);

		return *this;
	}


	// Get a command line with the program and its arguments, like you'd type into a shell or pass to CreateProcess()
	// Arguments with spaces will be double quoted.
	std::string getCommandlineString()
	{
		std::stringstream cmdline;

		cmdline << quoteIfNecessary(_program);

		for(std::vector<std::string>::iterator arg = _arguments.begin(); arg != _arguments.end(); ++arg)
		{
			cmdline << " " << quoteIfNecessary(*arg);
		}

		return cmdline.str();
	}

	// Execute the command and arguments, setting its standard streams to the three provided if they are not 0.
	// Blocks until it finishes.
	// If the command is not an absolute path, it will be searched for in the system PATH.
	// Returns the exit code of the process, or -1 if the process could not be started.

	int executeAndWait(int sin, int sout, int serr)
	{
	#ifdef WIN32
		STARTUPINFO startInfo;

		ZeroMemory(&startInfo, sizeof(startInfo));
		startInfo.cb = sizeof(startInfo);
		startInfo.dwFlags = STARTF_USESTDHANDLES;

		// convert file descriptors to win32 handles
		if(sin != 0)
		{
			startInfo.hStdInput = (HANDLE)_get_osfhandle(sin);
		}
		if(sout != 0)
		{
			startInfo.hStdOutput = (HANDLE)_get_osfhandle(sout);
		}
		if(serr != 0)
		{
			startInfo.hStdError = (HANDLE)_get_osfhandle(serr);
		}

		PROCESS_INFORMATION procInfo;
		if(CreateProcessA(NULL, const_cast<char*>(getCommandlineString().c_str()), NULL, NULL, true, 0, NULL, NULL, &startInfo, &procInfo) == 0)
		{
			int lasterror = GetLastError();

			LPTSTR strErrorMessage = NULL;

			FormatMessage(
				FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS | FORMAT_MESSAGE_ARGUMENT_ARRAY | FORMAT_MESSAGE_ALLOCATE_BUFFER,
				NULL,
				lasterror,
				0,
				(LPTSTR)(&strErrorMessage),
				0,
				NULL);

			std::cerr << "CreateProcess(" << getCommandlineString() << ") failed with error " << std::dec << lasterror << ": " << strErrorMessage << std::endl;
			return -1;
		}

		// Wait until child process exits.
		WaitForSingleObject(procInfo.hProcess, INFINITE);

		DWORD exitCode;
		GetExitCodeProcess(procInfo.hProcess, &exitCode);

		// Close process and thread handles.
		CloseHandle(procInfo.hProcess);
		CloseHandle(procInfo.hThread);

		return exitCode;
#else
	/*                   from David A Curry,
						 "Using C on the UNIX System" p. 105-106    */

		//create array of C strings
		std::vector<char*> execvpArguments;

		char* program_c_string = new char[_program.size() + 1];
		strcpy(program_c_string, _program.c_str());
		execvpArguments.push_back(program_c_string);

		for(std::vector<std::string>::iterator arg = _arguments.begin(); arg != _arguments.end(); ++arg)
		{
			char* c_string = new char[(*arg).size() + 1];
			strcpy(c_string, arg->c_str());
			execvpArguments.push_back(c_string);
		}

		//null-terminate the argument array
		execvpArguments.push_back(NULL);

		//we use pipes as the delimiter between logical arguments
//		std::cerr << "Invoking command: ";
//		for(char* arg : execvpArguments)
//		{
//			std::cerr << '|' << arg;
//		}
//		std::cerr << std::endl;

		int status;
		pid_t pid;

		if ((pid = fork()) < 0) {
			perror("fork");
			return -1;
		}

		if(pid == 0) {
			if( sin != 0 ) {
				close( 0 );  dup( sin );
			}
			if( sout != 1 ) {
				close( 1 );  dup( sout );
			}
			if( serr != 2 ) {
				close( 2 );  dup( serr );
			}

			execvp(_program.c_str(), &execvpArguments[0]);
			perror(("Error executing " + _program).c_str());
			exit(1);
		}

		//free arg strings
		for(std::vector<char*>::iterator arg = execvpArguments.begin(); arg != execvpArguments.end(); ++arg)
		{
			delete[] *arg;
		}

		while( wait(&status) != pid ) ;
		return( status );
	#endif
	}
};


int main( int argc, char *argv[] )
{
	
	//first, sort out amberhome
	if(getenv("AMBERHOME") == NULL )
	{
	   std::cerr << "error: AMBERHOME is not set!" << std::endl;
	   exit(1);
	}
	amberhome = getenv("AMBERHOME");
	
	std::string ofname;
	std::vector<std::string> cppArgs;
	
	//CC and LD are defined at build-time to the name of the compiler and linker
	std::string c_compiler(CC);
	std::string linker(LD);
	
	//named of input files without any file extension
	std::vector<std::string> inputFileBasenames;

	std::vector<std::string> inputObjectFiles;

	// parse nab options:
	ofname = "";
	for(int ac = 1; ac < argc; ac++ )
	{
		std::string currentArg = argv[ac];
		
		if(currentArg == "-c")
		{
			compileOnly = true;
		}
		else if(currentArg == "-cgd")
		{
			verbose = true;
			nab2cDebugMode = true;
		}
		else if(currentArg.substr(0, 2) == "-D" || currentArg.substr(0, 2) == "-I")
		{
			cppArgs.push_back(currentArg);
		}
		else if(currentArg == "--noassert" || currentArg == "-noassert")
		{
			noassert = true;
		}
		else if(currentArg == "--nodebug" || currentArg == "-nodebug")
		{
			nodebug = true;
		}
		else if(currentArg == "--save-intermediates")
		{
			saveIntermediates = true;
		}
		else if(currentArg == "-o")
		{
			ac++;
			if( ac == argc )
			{
				std::cerr << argv[ 0 ] << ": -o requires file name" << std::endl;
				return 1;
			}
			else
			{
				ofname = argv[ac];
			}
		}
		else if(currentArg == "--compiler")
		{
			ac++;
			if( ac == argc )
			{
				std::cerr << argv[ 0 ] << ": --compiler requires command for, or path to, compiler" << std::endl;
				return 1;
			}
			else
			{
				c_compiler = argv[ac];
			}
		}
		else if(currentArg == "--linker")
				{
					ac++;
					if( ac == argc )
					{
						std::cerr << argv[ 0 ] << ": --linker requires command for, or path to, linker" << std::endl;
						return 1;
					}
					else
					{
						linker = argv[ac];
					}
				}
		else if(currentArg == "-v")
		{
			verbose = true;
		}
		else if(currentArg.substr(0, 1) == "-")
		{
			std::cerr << "Invalid option " << currentArg << std::endl;
			return 1;
		}
		else
		{
			if(currentArg.substr(currentArg.size() - 4, 5) == ".nab")
			{
				inputFileBasenames.push_back(currentArg.substr(0, currentArg.size() - 4));
			}
			else if(currentArg.substr(currentArg.size() - 4, 5) == ".obj" || currentArg.substr(currentArg.size() - 2, 5) == ".o")
			{
				inputObjectFiles.push_back(currentArg);
			}
			else
			{
				std::cerr << "Input file " << currentArg << " does not have a .nab or object file extension! It will be ignored." << std::endl;
			}
		}
	}
	
	if(verbose)
	{
		std::cout << "nab++ wrapper version 1.3" << std::endl;
	}
	
	if(argc == 1 || (inputFileBasenames.empty() && inputObjectFiles.empty()))
	{
		std::cerr << "usage: " << argv[0] << " [-c] [-cgd] [-Dcompiledef] [--noassert] [--nodebug] [--save-intermediates] [--compiler alternate-c-compiler] [--linker alternate-fortran-linker] [-o executable-name] [-v] nab file(s)" << std::endl;
		return 1;
	}
	
	if(ofname == "")
	{
		if(inputFileBasenames.size() > 0)
		{
			//default to name of first .nab file supplied
			ofname = inputFileBasenames.at(0);
		}
		else
		{
			//default GCC program name
			ofname = "a.out";
		}

	}

	// compile .nab files to .o files, if there are any
	if(inputFileBasenames.size() > 0)
	{
		compileNab(inputFileBasenames, cppArgs);

		compileC(c_compiler, inputFileBasenames);
	}

	if(!compileOnly)
	{
		std::vector<std::string> allObjectFiles;

		//add object files for .nab files we've compiled
		for(std::vector<std::string>::iterator inputFileBasename = inputFileBasenames.begin(); inputFileBasename != inputFileBasenames.end(); ++inputFileBasename)
		{
			allObjectFiles.push_back(*inputFileBasename + CREATED_OBJECT_SUFFIX);
		}

		//add object files provided on the command line
		allObjectFiles.insert(allObjectFiles.end(), inputObjectFiles.begin(), inputObjectFiles.end());

		linkC(linker, allObjectFiles, ofname);
	}

	return 0;
}


// Creates and opens a temporary file with the provided prefix
// Returns its file descriptor and path.
static std::pair<int, std::string> openTempFile(std::string prefix)
{
	int filedesc;
	std::string filepath;

#ifdef WIN32	
	static char tempFolderPath[MAX_PATH];
	char filePathCString[MAX_PATH];
	
	GetTempPath(MAX_PATH, tempFolderPath);
	GetTempFileName(tempFolderPath, prefix.c_str(), 0, (LPSTR)(&filePathCString));
	
	if(verbose)
	{
		std::cout << "Creating temporary file for preprocessed intemediate: " << filePathCString << std::endl;
	}
	
	//allow the processes we spawn to inherit our handles
	SECURITY_ATTRIBUTES securityAttributes;
	ZeroMemory(&securityAttributes, sizeof(securityAttributes));

	securityAttributes.nLength = sizeof(SECURITY_ATTRIBUTES);
	securityAttributes.bInheritHandle = true;
	
	HANDLE fileHandle = CreateFileA(
		filePathCString,
		(GENERIC_READ | GENERIC_WRITE),
		FILE_SHARE_WRITE | FILE_SHARE_READ,
		&securityAttributes,
		CREATE_ALWAYS,
		FILE_ATTRIBUTE_NORMAL,
		NULL);
	
    if(fileHandle == INVALID_HANDLE_VALUE){
		std::cerr << "can't create temporary file " << filePathCString << "error: " << strerror(errno) << std::endl;
		exit(1);
	}
	
	filedesc = _open_osfhandle((intptr_t)fileHandle, _O_TEXT);
	filepath = filePathCString;
	
#else 
	std::string tempFilePathTemplate = "/tmp/" + prefix + "_XXXXXX";
	
	// needs to be mutable for mkstemp, so we copy the data out of the string
	char* filepath_cstring = (char*) malloc(tempFilePathTemplate.size() + 1);
	strcpy(filepath_cstring, tempFilePathTemplate.c_str());
	
	if( ( filedesc = mkstemp( filepath_cstring ) ) < 0 ){
		std::cerr << "can't create temporary file " << filepath_cstring << "error: " << strerror(errno) << std::endl;
		exit(1);
	}
	
	filepath = filepath_cstring;
	free(filepath_cstring);
#endif 

	return std::pair<int, std::string>(filedesc, filepath);
}

static void compileNab(std::vector<std::string> inputFileBasenames, std::vector<std::string> cppArgs)
{
	int	status;

	for(std::vector<std::string>::iterator inputFileBasename = inputFileBasenames.begin(); inputFileBasename != inputFileBasenames.end(); ++inputFileBasename)
	{

		//for this section we want to use the .nab extensions
		std::string inputFile = *inputFileBasename + ".nab";
		
		if(access(inputFile.c_str(), F_OK ) != 0)
		{
			std::cerr << "error: " << inputFile << ": no such file." << std::endl;
			exit(1);
		}

		//get temp file for cpp output, and run CPP
		//CPP is defined at compile time to be the preprocessor executable name
		CommandLine cppCmd(CPP);
		cppCmd.arg("-l");

		for(std::vector<std::string>::iterator cpparg = cppArgs.begin(); cpparg != cppArgs.end(); ++cpparg)
		{
			cppCmd.arg(*cpparg);
		}

		cppCmd.arg("-I").arg(amberhome + "/include").arg(inputFile);
		
		if(verbose)
		{
			std::cout << "Invoking preprocessor: " << cppCmd.getCommandlineString() << std::endl;
		}
		
		std::pair<int, std::string> tempFileData = openTempFile("cpp_of");
		
		int preprocessedTempFd = tempFileData.first;
		
		status = cppCmd.executeAndWait(0, preprocessedTempFd, STDERR_FILENO);
		
		if( status != 0)
		{
			if(status == -1)
			{
				std::cerr << "Failed to invoke the preprocessor." << std::endl;
			}
			else
			{
				std::cerr << "Preprocessor failed with exit code " << status << std::endl;
			}
			exit(2);
		}
		
		lseek(preprocessedTempFd, 0L, SEEK_SET);

		// next run nab2c 
#ifdef MPI
		CommandLine nabCmd(amberhome + "/bin/mpinab2c");
#else
		CommandLine nabCmd(amberhome + "/bin/nab2c");
#endif
	
		if(nab2cDebugMode)
		{
			nabCmd.arg("-cgdebug");
		}
		if(noassert)
		{
			nabCmd.arg("-noassert");
		}
		if(nodebug)
		{
			nabCmd.arg("-nodebug");
		}
		
		nabCmd.arg("-nfname").arg(inputFile);
		
		if(verbose)
		{
			std::cout << "Invoking C generator: " << nabCmd.getCommandlineString() << std::endl;
		}
		
		//this generates a .c file in the current directory corresponding to the .nab input file
		status = nabCmd.executeAndWait(preprocessedTempFd, STDOUT_FILENO, STDERR_FILENO);
		
		if( status != 0)
		{
			if(status == -1)
			{
				std::cerr << "Failed to invoke the C generator." << std::endl;
			}
			else
			{
				std::cerr << "C generator failed with exit code " << status << std::endl;
			}
			exit(2);
		}
		
		close(preprocessedTempFd);
		
		//delete temp file
		if(!saveIntermediates)
		{
			if(remove(tempFileData.second.c_str()) != 0)
			{
				perror("remove");
			}
		}
	}
}

// compile all if the C files generated by compileNab into .o files
static void compileC(std::string c_compiler, std::vector<std::string> inputFileBasenames)
{
	for(std::vector<std::string>::iterator inputFileBasename = inputFileBasenames.begin(); inputFileBasename != inputFileBasenames.end(); ++inputFileBasename)
	{
		CommandLine cmd(c_compiler);

		cmd.arg("-I" + amberhome + "/include");

		cmd.arg("-c");

		cmd.arg(*inputFileBasename + ".c");


		for(std::vector<std::string>::iterator flag = builtinCompilerFlags.begin(); flag != builtinCompilerFlags.end(); ++flag)
		{
			cmd.arg(*flag);
		}

		cmd.arg("-o").arg(*inputFileBasename + CREATED_OBJECT_SUFFIX);

		if(verbose)
		{
			std::cout << "Invoking compiler: " << cmd.getCommandlineString() << std::endl;
		}

		int status = cmd.executeAndWait(STDIN_FILENO, STDOUT_FILENO, STDERR_FILENO);
		if( status != 0)
		{
			if(status == -1)
			{
				std::cerr << "Failed to invoke the C compiler." << std::endl;
			}
			else
			{
				std::cerr << "C compiler failed with exit code " << status << std::endl;
			}
			exit(2);
		}
	}
	
	//delete temp file
	if(!saveIntermediates)
	{
		for(std::vector<std::string>::iterator inputFileBasename = inputFileBasenames.begin(); inputFileBasename != inputFileBasenames.end(); ++inputFileBasename)
		{
			if(remove((*inputFileBasename + ".c").c_str()) != 0)
			{
				perror("remove");
			}		
		}
	}
}

// Link all of the object files provided into an executable
static void linkC(std::string linker, std::vector<std::string> objectFiles, std::string outputFilename)
{
	CommandLine cmd(linker);

	for(std::vector<std::string>::iterator objectFile = objectFiles.begin(); objectFile != objectFiles.end(); ++objectFile)
	{
		cmd.arg(*objectFile);
	}

	for(std::vector<std::string>::iterator flag = builtinCompilerFlags.begin(); flag != builtinCompilerFlags.end(); ++flag)
	{
		cmd.arg(*flag);
	}

	cmd.arg("-L" + amberhome + "/lib");

	// set up RPATH
#if defined(__APPLE__) || defined(__linux__)
	// On OSX, Amber libraries use '@rpath/<name>.dylib` as their install_name, so we have to add $AMBERHOME/lib as the RPATH
	// On Linux, we can still fall back to LD_LIBRARY_PATH, but setting the RPATH makes it a bit more convenient for the user.
	cmd.arg("-Wl,-rpath," + amberhome + "/lib");
#endif

	for(std::vector<std::string>::iterator library = builtinLibraries.begin(); library != builtinLibraries.end(); ++library)
	{

		cmd.arg(*library);
	}

	for(std::vector<std::string>::iterator linkDir = builtinLinkDirs.begin(); linkDir != builtinLinkDirs.end(); ++linkDir)
	{
		cmd.arg("-L" + *linkDir);

		// also add RPATHs on Linux
#if defined(__linux__)
		cmd.arg("-Wl,-rpath," + *linkDir);
#endif

	}

	cmd.arg("-o").arg(outputFilename);

	if(verbose)
	{
		std::cout << "Invoking linker: " << cmd.getCommandlineString() << std::endl;
	}

	int8_t status = cmd.executeAndWait(STDIN_FILENO, STDOUT_FILENO, STDERR_FILENO);
	if( status != 0)
	{
		if(status == -1)
		{
			std::cerr << "Failed to invoke the linker." << std::endl;
		}
		else
		{
			std::cerr << "linker failed with exit code " << status << std::endl;
		}
		exit(2);
	}

	//delete temp file
	if(!saveIntermediates)
	{
		for(std::vector<std::string>::iterator objectFile = objectFiles.begin(); objectFile != objectFiles.end(); ++objectFile)
		{
			if(remove(objectFile->c_str()) != 0)
			{
				perror("remove");
			}
		}
	}
}
