#include <ql/qldefines.hpp>
#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#    include <ql/auto_link.hpp>
#endif
#include <ql/currencies/europe.hpp>
#include <ql/experimental/credit/cdo.hpp>
#include <ql/experimental/credit/gaussianlhplossmodel.hpp>
#include <ql/experimental/credit/homogeneouspooldef.hpp>
#include <ql/experimental/credit/inhomogeneouspooldef.hpp>
#include <ql/experimental/credit/integralcdoengine.hpp>
#include <ql/experimental/credit/midpointcdoengine.hpp>
#include <ql/experimental/credit/pool.hpp>
#include <ql/experimental/credit/randomdefaultlatentmodel.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/credit/flathazardrate.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <iomanip>
#include <iostream>

// This makes it easier to use array literals (for new code, use std::vector though)
#define LENGTH(a) (sizeof(a)/sizeof(a[0]))