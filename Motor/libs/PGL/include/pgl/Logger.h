#pragma once

#ifdef Mt_Debug
#pragma message("Mt_Debug is defined")
#endif

#ifdef Mt_Release
#pragma message("Mt_Release is defined")
#endif


#include <iostream>
#include <ostream>
#include <ctime>
#include <string>


struct Logger_Color
{
	enum class Colors
	{
		RED,
		GREEN,
		BLUE,
		YELLOW,
		WHITE,
		BLACK,
		PURPLE,
		BROWN
	};

	Colors font_color;
	Colors background_color;

	Logger_Color(Colors fColor = Colors::WHITE, Colors bckColor = Colors::BLACK)
		: font_color(fColor), background_color(bckColor)
	{}

};


class Logger
{
public:

	enum Log_Level : char
	{
        CRIT 			= 0,
		ERROR			= 1,
		WARN			= 2,
		INFO			= 3
	};

	enum Filter : char
	{
		FILTER_CRIT 		    = 0,
		FILTER_ERROR 	 		= 1,
		FILTER_WARN				= 2,
		FILTER_INFO				= 3
	};



	Logger(Filter filter, Logger_Color _color = Logger_Color{ Logger_Color::Colors::WHITE, Logger_Color::Colors::BLACK});
	Logger(const char* filepath);

	~Logger();
	
	void Log(const char* msg, Log_Level lvl = INFO, Logger_Color color = Logger_Color{ Logger_Color::Colors::WHITE, Logger_Color::Colors::BLACK }) const;
	
	void SetColor(Logger_Color color) const;
	void SetFilterLevel(Filter filter);
	void SetFilepath(const char* filepath);


	const char* Filepath() const { return m_Filepath; }
	Filter filter() const { return m_Filter; }
	Logger_Color color() const { return m_Color; }

private :

	Logger_Color m_Color;

	bool m_InTerminal;

	bool m_InFile;
	const char* m_Filepath;

	Filter m_Filter;
};


// inline functions must be in the same translation unit
inline void Logger::Log(const char* msg, Logger::Log_Level lvl, Logger_Color color) const
{
	#ifdef Mt_Debug

	unsigned char _lvl = lvl;
	unsigned char _filter = m_Filter;

	

	auto time = std::time(nullptr);
	//auto date = std::asctime(std::localtime(&time));

	char date[26];
	ctime_s(date, sizeof date, &time);
	
	std::string prefix = std::string("\n") + date + std::string("[ line : ") + std::to_string(__LINE__) + "]" + std::string(", [file : ") + __FILE__ + std::string("]     message :\n");

	// we cut out messages at two high level
	if ((char)lvl > (char)m_Filter)
	{

		Logger::SetColor(Logger_Color::Colors::YELLOW);

		std::cout << prefix << "cut message : too low level\n";

		// back to default
		std::cout << "\033[0m";
		Logger::SetColor(m_Color);

		return;
	}


	switch (lvl)
	{
	case Logger::CRIT:
		// set color
		Logger::SetColor(color);

		std::cerr << "\n\n***************   CRITICAL   ***************\n";

		// print msg	 
		std::cerr << prefix << msg;

		std::cerr << "\n********************************************\n\n";

		// back to default
		std::cerr << "\033[0m";
		Logger::SetColor(m_Color);
		break;


	case Logger::ERROR:
		// set color
		Logger::SetColor(color);

		// print msg
		std::cerr << prefix << msg;

		// back to default
		std::cerr << "\033[0m";
		Logger::SetColor(m_Color);

		break;

	case Logger::WARN:
		// set color
		Logger::SetColor(color);

		// print msg
		std::cout << prefix << msg;

		// back to default
		std::cout << "\033[0m";
		Logger::SetColor(m_Color);
		break;

	case Logger::INFO:
		// set color
		Logger::SetColor(color);

		// print msg
		std::cout << prefix << msg;

		// back to default
		std::cout << "\033[0m";
		Logger::SetColor(m_Color);
		break;

	default:
		std::cout << prefix << "Errors in switch in Log function\n";
		break;
	}
	#else // Mt_Debug

	return;

	#endif // Mt_Debug	
}