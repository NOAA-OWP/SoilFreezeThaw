#ifndef HPP_STRING_VECBUF
#define HPP_STRING_VECBUF

#include <streambuf>
#include <vector>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>

using vecbuf = std::vector<char>;
using OStreamType = boost::iostreams::stream<boost::iostreams::back_insert_device<vecbuf>>;

class membuf : public std::streambuf {
public:
    membuf(char *begin, size_t size) {
        this->setg(begin, begin, begin + size);
    }
};

#endif
