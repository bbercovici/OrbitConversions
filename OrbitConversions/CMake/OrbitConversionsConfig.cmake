# MIT License
# Copyright (c) 2017 Benjamin Bercovici

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
if (EXISTS /home/bebe0705/.am_fortuna)
	set(OC_INCLUDE_DIR /home/bebe0705/libs/local/include/)
	set(OC_LIBRARY /home/bebe0705/libs/local/lib/libOrbitConversions.so)
else()

	if (APPLE)
		set(OC_LIBRARY /usr/local/lib/libOrbitConversions.dylib)
	elseif(UNIX AND NOT APPLE)
		set(OC_LIBRARY /usr/local/lib/libOrbitConversions.so)
	else()
		message(FATAL_ERROR "Unsupported platform")
	endif()
	set(OC_INCLUDE_DIR /usr/local/include/OrbitConversions/)
endif()

message("-- Found OrbitConversions: " ${OC_LIBRARY})
