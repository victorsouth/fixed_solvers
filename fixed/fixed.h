#pragma once

#ifndef __FIXED_H__
#define __FIXED_H__

#define _USE_MATH_DEFINES // подключить константы
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fixed/helpers/math_helpers.h"

namespace fixed_solvers {
;


inline std::string wide2string(const std::wstring& str)
{
    std::string ascii_string(str.size(), '\0');
    std::transform(str.begin(), str.end(), ascii_string.begin(),
        [](wchar_t s) { return static_cast<char>(s); });
    return ascii_string;
}

/// from https://stackoverflow.com/questions/148403/utf8-to-from-wide-char-conversion-in-stl/148665#148665
inline std::wstring UTF8_to_wchar(const char* in)
{
    std::wstring out;
    unsigned int codepoint;
    while (*in != 0)
    {
        unsigned char ch = static_cast<unsigned char>(*in);
        if (ch <= 0x7f)
            codepoint = ch;
        else if (ch <= 0xbf)
            codepoint = (ch << 6) | (ch & 0x3f);
        else if (ch <= 0xdf)
            codepoint = ch & 0x1f;
        else if (ch <= 0xef)
            codepoint = ch & 0x0f;
        else
            codepoint = ch & 0x07;
        ++in;
        if (((*in & 0xc0) != 0x80) && (codepoint <= 0x10ffff))
        {
            if constexpr (sizeof(wchar_t) > 2)
                out.append(1, static_cast<wchar_t>(codepoint));
            else if (codepoint > 0xffff)
            {
                out.append(1, static_cast<wchar_t>(0xd800 + (codepoint >> 10)));
                out.append(1, static_cast<wchar_t>(0xdc00 + (codepoint & 0x03ff)));
            }
            else if (codepoint < 0xd800 || codepoint >= 0xe000)
                out.append(1, static_cast<wchar_t>(codepoint));
        }
    }
    return out;
}

inline std::wstring string2wide(const std::string& str)
{
#ifndef __GNUC__
    try {
        std::wstring wstr(str.size(), 0);
        std::use_facet<std::ctype<wchar_t> >(std::locale("rus_rus.1251")).widen
        (&str[0], &str[0] + str.size(), &wstr[0]);
        return wstr;
    }
    catch (...)
    {
        std::wstring wide_string(str.begin(), str.end());
        return wide_string;
    }
#else
    std::wstring wstr = UTF8_to_wchar(str.c_str());
    return wstr;
#endif
}


template<class T>
inline std::string int2str(T i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

template<class T>
inline std::wstring int2wstr(T i)
{
    std::wstringstream ss;
    ss << i;
    return ss.str();
}

inline void string_replace(std::string& str, const std::string& from, const std::string& to)
{
    if (from.empty())
        return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}


}

#include <iomanip>
#include "fixed/qp/qp_wrapper.h"
#include "fixed/array_ext.h"
#include "fixed/fixed_system.h"
#include "fixed/fixed_linear_solver.h"
#include "fixed/fixed_constraints.h"
#include "fixed/fixed_nonlinear_solver.h"
#include "fixed/fixed_optimizer.h"

#endif
