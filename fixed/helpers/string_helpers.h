#pragma once

#include <string>
#include <locale>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <boost/type_index.hpp>

namespace fixed_solvers {

inline bool string_ends_with(std::string const& value, std::string const& ending)
{
    if (ending.size() > value.size())
        return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

/// @brief from https://stackoverflow.com/questions/148403/utf8-to-from-wide-char-conversion-in-stl/148665#148665
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


/// @brief Преобразует стандартную строку (std::string) в широкую строку (std::wstring).
/// Функция корректно обрабатывает строки в кодировке UTF-8 или системной локали.
/// На Windows используется локаль "rus_rus.1251" для расширения символов.
/// В GCC/Linux возвращается прямое преобразование через UTF-8.
/// @param str Входная строка типа std::string.
/// @return std::wstring — преобразованная строка с символами wchar_t.
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


inline std::string wide2string(const std::wstring& str)
{
    std::string ascii_string(str.size(), '\0');
    std::transform(str.begin(), str.end(), ascii_string.begin(),
        [](wchar_t s) { return static_cast<char>(s); });
    return ascii_string;
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

/// @brief Получает имя класса из типа объекта
/// @param t объект для определения типа
/// @return строка с именем класса
#ifdef _MSC_VER
template<typename T>
inline std::string __CLASS__(const T& t) {
    std::string str = boost::typeindex::type_id<decltype(t)>().pretty_name();
    auto ppos = str.find_last_of(" \t");
    if (ppos != str.npos)str = str.substr(ppos + 1);
    return str;
}
#else
#undef __CLASS__  // Убираем макрос, если он был определен ранее
template<typename T>
inline std::string __CLASS__(const T& t) {
    std::string str = boost::typeindex::type_id<decltype(t)>().pretty_name();
    return str;
}
#endif

/// @brief Сохраняет вектор в поток
/// @param os выходной поток
/// @param v вектор для сохранения
/// @return ссылка на поток
template <typename T>
inline std::ostream& save_vector(std::ostream& os, const std::vector<T>& v) {
    os << v.size() << '\n';
    for (auto val : v) {
        os << val << '\n';
    }
    return os << '\n';
}

/// @brief Загружает вектор из потока
/// @param is входной поток
/// @param v вектор для загрузки
/// @return ссылка на поток
template <typename T>
inline std::istream& load_vector(std::istream& is, std::vector<T>& v) {
    size_t v_size;
    is >> v_size;
    if (is.fail()) {
        throw std::runtime_error("bad v_size");
    }
    v.resize(v_size);
    for (auto& val : v) {
        is >> val;
    }
    if (is.fail()) {
        throw std::runtime_error("v read error");
    }
    return is;
}


}

