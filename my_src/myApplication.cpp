/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Jose Aparicio
 Copyright (C) 2014 Peter Caspers

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include "myApplication.hpp"

using namespace std;
using namespace QuantLib;


namespace my_application_cdo_test {

    // 4 tranches
    // 1st tranche: 0-3%
    // 2nd tranche: 3-6%
    // 3rd tranche: 6-10%
    // 4th tranche: 10-100%
    // nth tranche: [hwAttachment[n-1], hwDetachment[n-1]]
    // Real hwAttachment[] = {0.00, 0.03, 0.06, 0.10};
    // Real hwDetachment[] = {0.03, 0.06, 0.10, 1.00};



    struct hwDatum {
        Real correlation; // correlation between the defaults of the reference entities in the
                          // Synthetic CDO. i.e. the likelihood that the entities will default on
                          // their debt obligations at the same time. Note that this is under
                          // "single-factor model", assumed that the default of each entity is
                          // driven by a single common factor, along with an idiosyncratic
                          // (entity-specific) factor.

        Integer
            nm; // degrees of freedom parameter for a Student's t-distribution in a t-copula model. if -1, then gaussian
        Integer
            nz; // degrees of freedom parameter for a Student's t-distribution in a t-copula model. if -1, then gaussian

        Real trancheSpread[4]; // expected tranche spreads for the 4 tranches. used in the test
                               // cases to compare the results of the pricing engines.
    };


    // HW Table 7
    // corr, Nm, Nz, 0-3, 3-6, 6-10, 10-100
    hwDatum hwData7[] = {{0.1, -1, -1, {2279, 450, 89, 1}},
                         {0.3, -1, -1, {1487, 472, 203, 7}},
                         // Opening the T, T&G tests too. The convolution is analytical
                         //   now so it runs it a time comparable to the gaussian tests and
                         //   has enough precission to pass the tests.
                         // Below the T models are integrated with a quadrature, even if this
                         //   is incorrect the test pass good enough, the quadrature gets to
                         //   be worst as the kernel deviates from a normal, this is low
                         //   orders of the T; here 5 is enough, 3 would not be.
                         {0.3, -1, 5, {1766, 420, 161, 6}},
                         {0.3, 5, -1, {1444, 408, 171, 10}},
                         {0.3, 5, 5, {1713, 359, 136, 9}}};


    struct tranches {
        Real attachment;
        Real detachment;
    };


    // 4 tranches
    // 1st tranche: 0-3%
    // 2nd tranche: 3-6%
    // 3rd tranche: 6-10%
    // 4th tranche: 10-100%
    // tranches hwTranches[] = {{0.00, 0.03}, {0.03, 0.06}, {0.06, 0.10}, {0.10, 1.00}};
    vector<tranches> hwTranches = {{0.00, 0.03}, {0.03, 0.06}, {0.06, 0.10}, {0.10, 1.00}};


    void check(int i,
               int j,
               const std::string& desc,
               Real found,
               Real expected,
               Real bpTolerance,
               Real relativeTolerance) {

        Real absDiff = found - expected;
        Real relDiff = absDiff / expected;
        cout << "case " << i << " " << j << " (" << desc << "): " << found << " vs. " << expected << endl;
    }


    void runTest(unsigned dataSet) {
        cout << endl;
        cout << "Testing CDO premiums against Hull-White values for data set " << dataSet << "..." << endl;
        using namespace my_application_cdo_test;

#pragma region parameters
        Size poolSize = 100; // number of entities in the basket

        Real lambda = 0.01; // hazard rate, i.e. the instantaneous risk of default at a given time
                            // t, given that no default has occurred up until that time

        // [model param]
        // nBuckets and period determine the computation time
        Size nBuckets = 200;

        // [model param]
        // Period period = 1*Months;
        // for MC engines
        Size numSims = 5000;

        // [for construction of the risk free yield term structure]
        Real rate = 0.05;

        // [for construction of the risk free yield term structure]
        DayCounter daycount = Actual360();

        // [for construction of the risk free yield term structure]
        Compounding cmp = Continuous; // Simple;

        Real recovery = 0.4;
        std::vector<Real> nominals(poolSize, 100.0);
        Real premium = 0.02;
        Period maxTerm(5, Years);
        Schedule schedule = MakeSchedule()
                                .from(Date(1, September, 2006))
                                .to(Date(1, September, 2011))
                                .withTenor(Period(3, Months))
                                .withCalendar(TARGET());

        Date asofDate = Date(31, August, 2006);
#pragma endregion

        Settings::instance().evaluationDate() = asofDate;

        //  risk-free rate term structure
        ext::shared_ptr<YieldTermStructure> yieldPtr(new FlatForward(asofDate, rate, daycount, cmp));
        Handle<YieldTermStructure> yieldHandle(yieldPtr);

        // hazard rate
        Handle<Quote> hazardRate(ext::shared_ptr<Quote>(new SimpleQuote(lambda)));

        std::vector<Handle<DefaultProbabilityTermStructure>> basket;
        ext::shared_ptr<DefaultProbabilityTermStructure> ptr(
            new FlatHazardRate(asofDate, hazardRate, ActualActual(ActualActual::ISDA)));
        ext::shared_ptr<Pool> pool(new Pool());
        std::vector<std::string> names;
        // probability key items
        std::vector<Issuer> issuers;
        std::vector<std::pair<DefaultProbKey, Handle<DefaultProbabilityTermStructure>>> probabilities;
        probabilities.emplace_back(
                                    NorthAmericaCorpDefaultKey(EURCurrency(), SeniorSec, Period(0, Weeks), 10.),
                                    Handle<DefaultProbabilityTermStructure>(ptr)
                                );

        for (Size i = 0; i < poolSize; ++i) {
            std::ostringstream o;
            o << "issuer-" << i;
            names.push_back(o.str());
            basket.emplace_back(ptr);
            issuers.emplace_back(probabilities);
            pool->add(names.back(), issuers.back(),
                      NorthAmericaCorpDefaultKey(EURCurrency(), QuantLib::SeniorSec, Period(), 1.));
        }

        ext::shared_ptr<SimpleQuote> correlation(new SimpleQuote(0.0));
        Handle<Quote> hCorrelation(correlation);
        // QL_REQUIRE(LENGTH(hwAttachment) == LENGTH(hwDetachment), "data length does not match");

        ext::shared_ptr<PricingEngine> midPCDOEngine(new MidPointCDOEngine(yieldHandle));
        ext::shared_ptr<PricingEngine> integralCDOEngine(new IntegralCDOEngine(yieldHandle));

        const Size i = dataSet;
        correlation->setValue(hwData7[i].correlation);
        // QL_REQUIRE(LENGTH(hwAttachment) == LENGTH(hwData7[i].trancheSpread),
        //            "data length does not match");


// define the models to be tested
// update the tolerance values based on the models
#pragma region Models
        std::vector<ext::shared_ptr<DefaultLossModel>> basketModels;
        std::vector<std::string> modelNames;
        std::vector<Real> relativeToleranceMidp, relativeTolerancePeriod, absoluteTolerance;

        if (hwData7[i].nm == -1 && hwData7[i].nz == -1) {
            ext::shared_ptr<GaussianConstantLossLM> gaussKtLossLM(
                new GaussianConstantLossLM(hCorrelation, std::vector<Real>(poolSize, recovery),
                                           LatentModelIntegrationType::GaussianQuadrature, poolSize,
                                           GaussianCopulaPolicy::initTraits()));

            // 1.-Inhomogeneous gaussian
            modelNames.emplace_back("Inhomogeneous gaussian");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new IHGaussPoolLossModel(gaussKtLossLM, nBuckets, 5., -5, 15)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.04);
            relativeTolerancePeriod.push_back(0.04);
            // 2.-homogeneous gaussian
            modelNames.emplace_back("Homogeneous gaussian");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new HomogGaussPoolLossModel(gaussKtLossLM, nBuckets, 5., -5, 15)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.04);
            relativeTolerancePeriod.push_back(0.04);
            // 3.-random default gaussian
            modelNames.emplace_back("Random default gaussian");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new RandomDefaultLM<GaussianCopulaPolicy>(gaussKtLossLM, numSims)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.07);
            relativeTolerancePeriod.push_back(0.07);
            // SECOND MC
            // gaussian LHP
            modelNames.emplace_back("Gaussian LHP");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new GaussianLHPLossModel(hCorrelation, std::vector<Real>(poolSize, recovery))));
            absoluteTolerance.push_back(10.);
            relativeToleranceMidp.push_back(0.5);
            relativeTolerancePeriod.push_back(0.5);
            // Binomial...
            // Saddle point...
            // Recursive ...
        } else if (hwData7[i].nm > 0 && hwData7[i].nz > 0) {
            TCopulaPolicy::initTraits initTG;
            initTG.tOrders.push_back(hwData7[i].nm);
            initTG.tOrders.push_back(hwData7[i].nz);
            ext::shared_ptr<TConstantLossLM> TKtLossLM(new TConstantLossLM(
                hCorrelation, std::vector<Real>(poolSize, recovery),
                LatentModelIntegrationType::GaussianQuadrature, poolSize, initTG));
            // 1.-inhomogeneous studentT
            modelNames.emplace_back("Inhomogeneous student");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new IHStudentPoolLossModel(TKtLossLM, nBuckets, 5., -5., 15)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.04);
            relativeTolerancePeriod.push_back(0.04);
            // 2.-homogeneous student T
            modelNames.emplace_back("Homogeneous student");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new HomogTPoolLossModel(TKtLossLM, nBuckets, 5., -5., 15)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.04);
            relativeTolerancePeriod.push_back(0.04);
            // 3.-random default student T
            modelNames.emplace_back("Random default studentT");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new RandomDefaultLM<TCopulaPolicy>(TKtLossLM, numSims)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.07);
            relativeTolerancePeriod.push_back(0.07);
            // SECOND MC
            // Binomial...
            // Saddle point...
            // Recursive ...
        } else if (hwData7[i].nm > 0 && hwData7[i].nz == -1) {
            TCopulaPolicy::initTraits initTG;
            initTG.tOrders.push_back(hwData7[i].nm);
            initTG.tOrders.push_back(45);
            /* T_{55} is pretty close to a gaussian. Probably theres no need to
            be this conservative as the polynomial convolution gets shorter and
            faster as the order decreases.
            */
            ext::shared_ptr<TConstantLossLM> TKtLossLM(new TConstantLossLM(
                hCorrelation, std::vector<Real>(poolSize, recovery),
                LatentModelIntegrationType::GaussianQuadrature, poolSize, initTG));
            // 1.-inhomogeneous
            modelNames.emplace_back("Inhomogeneous student-gaussian");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new IHStudentPoolLossModel(TKtLossLM, nBuckets, 5., -5., 15)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.04);
            relativeTolerancePeriod.push_back(0.04);
            // 2.-homogeneous
            modelNames.emplace_back("Homogeneous student-gaussian");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new HomogTPoolLossModel(TKtLossLM, nBuckets, 5., -5., 15)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.04);
            relativeTolerancePeriod.push_back(0.04);
            // 3.-random default
            modelNames.emplace_back("Random default student-gaussian");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new RandomDefaultLM<TCopulaPolicy>(TKtLossLM, numSims)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.07);
            relativeTolerancePeriod.push_back(0.07);
            // SECOND MC
            // Binomial...
            // Saddle point...
            // Recursive ...
        } else if (hwData7[i].nm == -1 && hwData7[i].nz > 0) {
            TCopulaPolicy::initTraits initTG;
            initTG.tOrders.push_back(45); // pretty close to gaussian
            initTG.tOrders.push_back(hwData7[i].nz);
            ext::shared_ptr<TConstantLossLM> TKtLossLM(new TConstantLossLM(
                hCorrelation, std::vector<Real>(poolSize, recovery),
                LatentModelIntegrationType::GaussianQuadrature, poolSize, initTG));
            // 1.-inhomogeneous gaussian
            modelNames.emplace_back("Inhomogeneous gaussian-student");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new IHStudentPoolLossModel(TKtLossLM, nBuckets, 5., -5., 15)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.04);
            relativeTolerancePeriod.push_back(0.04);
            // 2.-homogeneous gaussian
            modelNames.emplace_back("Homogeneous gaussian-student");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new HomogTPoolLossModel(TKtLossLM, nBuckets, 5., -5., 15)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.04);
            relativeTolerancePeriod.push_back(0.04);
            // 3.-random default gaussian
            modelNames.emplace_back("Random default gaussian-student");
            basketModels.push_back(ext::shared_ptr<DefaultLossModel>(
                new RandomDefaultLM<TCopulaPolicy>(TKtLossLM, numSims)));
            absoluteTolerance.push_back(1.);
            relativeToleranceMidp.push_back(0.07);
            relativeTolerancePeriod.push_back(0.07);
            // SECOND MC
            // Binomial...
            // Saddle point...
            // Recursive ...
        } else {
            return;
        }
#pragma endregion

        // start the valuation for each tranche
        for (Size j = 0; j < hwTranches.size(); j++) {

            ext::shared_ptr<Basket> basketPtr(new Basket(
                asofDate, // as of date, date, reference date for the basket. It could represent the
                          // date on which the basket was created or the date from which the credit
                          // risk of the entities in the basket is being measured.

                names, // names, list of strings, the names of the reference entities


                nominals, // nominals, list of real numbers, the notional amounts associated with
                          // each reference entity. The notional amount is the hypothetical amount
                          // of debt upon which the credit derivative's payments are based. It
                          // represents the face value of the debt associated with each reference
                          // entity.

                pool, // pool, share ptr, collection of reference entities such as corporate bonds
                      // or loans that the CDO is based on

                hwTranches[j].attachment, // attachment point, real number, the
                                          // threshold of losses that must be reached before the
                                          // protection seller (or the investors in the case of a
                                          // Synthetic CDO) starts to suffer losses.

                hwTranches[j].detachment    // detachment point, real number, the threshold of losses that must
                                            // be reached before the protection seller (or the investors in the
                                            // case of a Synthetic CDO) starts to suffer losses.
            ));

            std::ostringstream trancheId;

            trancheId << "Tranche #" << j + 1 << " "
                      << "[" << hwTranches[j].attachment << " - " << hwTranches[j].detachment << "]";

            SyntheticCDO cdoe(
                basketPtr, // basket, share ptr, collection of reference entities such as corporate bonds or
                           // loans that the CDO is based on

                Protection::Seller, // protection side, enum, the side of the credit default swap (CDS)
                                    // that the protection seller is on

                schedule,           // schedule, object, the schedule of payment dates for the CDO

                0.0,     // upfront payment rate, real number, a percentage of the notional amount that the
                         // protection buyer pays the protection seller at the start of the contract
                
                premium, // regular payment rate, real number, a percentage of the notional amount that the
                         // protection buyer pays to the protection seller on each payment date in
                         // the schedule
                
                daycount, // day counter,  used for calculating the time between dates
                
                Following // payment convention, convention for adjusting payment dates if they fall
                          // on a non-business day
            );

            for (Size im = 0; im < basketModels.size(); im++) {

                basketPtr->setLossModel(basketModels[im]);

                cdoe.setPricingEngine(midPCDOEngine);
                check(
                        i, // data set
                        j, // tranche
                        modelNames[im] + std::string(" with midp integration on ") + trancheId.str(), // description
                        cdoe.fairPremium() * 1e4, // found 
                        hwData7[i].trancheSpread[j], // expected
                        absoluteTolerance[im], // absolute tolerance
                        relativeToleranceMidp[im] // relative tolerance
                    );

                cdoe.setPricingEngine(integralCDOEngine);
                check(
                        i, // data set
                        j, // tranche
                        modelNames[im] + std::string(" with step integration on ") + trancheId.str(), // description
                        cdoe.fairPremium() * 1e4, // found
                        hwData7[i].trancheSpread[j], // expected
                        absoluteTolerance[im], // absolute tolerance
                        relativeTolerancePeriod[im] // relative tolerance
                    );
            }
        }
    }
}

int main(int argc, char* argv[]) {
    std::cout << "welcome to my application" << std::endl;
    my_application_cdo_test::runTest(0);
    my_application_cdo_test::runTest(1);
    my_application_cdo_test::runTest(2);
    my_application_cdo_test::runTest(3);
    my_application_cdo_test::runTest(4);
}
