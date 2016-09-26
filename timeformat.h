/*
	timeformat.h : This is a header FILE for the program of hash table.
    Copyright (C) 2016 by SHINEUKE ONO, JPN (given sur, country) <shinx55@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define MINUTE_SEC 60
#define HOUR_SEC (MINUTE_SEC * 60)//3,600
#define DAY_SEC (HOUR_SEC * 24)//86,400
#define WEEK_SEC (DAY_SEC * 7)//604,800
#define YEAR_SEC ((DAY_SEC * 365) + (DAY_SEC / 4)) //31,557,600
#define MONTH_SEC (YEAR_SEC / 12) //2,629,800

extern const char * formatSecond(char * a_buff, size_t a_buffSize, double a_timeSec, int a_pretty);
extern char * formatTime(char * a_timetxt, int a_sec);
extern int readTime(const char * a_buff);