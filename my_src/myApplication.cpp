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


    struct model{
        string modelName;
        ext::shared_ptr<DefaultLossModel> basketModel;
        Real absoluteTolerance;
        Real relativeToleranceMidp;
        Real relativeTolerancePeriod;


        // constructor
        model(const string& modelName,
          const ext::shared_ptr<DefaultLossModel>& basketModel,
          const Real& absoluteTolerance,
          const Real& relativeToleranceMidp,
          const Real& relativeTolerancePeriod)
        : modelName(modelName),
          basketModel(basketModel),
          absoluteTolerance(absoluteTolerance),
          relativeToleranceMidp(relativeToleranceMidp),
          relativeTolerancePeriod(relativeTolerancePeriod)
        {}
    };


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
        // used for random default model
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

        // basket
        std::vector<Handle<DefaultProbabilityTermStructure>> basket;

        ext::shared_ptr<DefaultProbabilityTermStructure> prob_term(
            new FlatHazardRate(asofDate, hazardRate, ActualActual(ActualActual::ISDA)));

        ext::shared_ptr<Pool> pool(new Pool());

        std::vector<std::string> names;

        // probability key items
        std::vector<Issuer> issuers;
        std::vector<std::pair<DefaultProbKey, Handle<DefaultProbabilityTermStructure>>> probabilities;
        probabilities.emplace_back(
                                    NorthAmericaCorpDefaultKey(EURCurrency(), SeniorSec, Period(0, Weeks), 10.),
                                    Handle<DefaultProbabilityTermStructure>(prob_term)
                                );

        for (Size i = 0; i < poolSize; ++i) {
            
            // construct the name of the reference entity and push it into the vector
            std::ostringstream o;
            o << "issuer-" << i;
            names.push_back(o.str());

            // add the probability term structure to the basket
            basket.emplace_back(prob_term);

            // add the probability key items to the vector
            issuers.emplace_back(probabilities);

            // add the reference entities to the pool
            pool->add(
                        names.back(), 
                        issuers.back(),
                        NorthAmericaCorpDefaultKey(EURCurrency(), QuantLib::SeniorSec, Period(), 1.)
                    );
        }

        ext::shared_ptr<SimpleQuote> correlation(new SimpleQuote(0.0));
        Handle<Quote> hCorrelation(correlation);

        ext::shared_ptr<PricingEngine> midPCDOEngine(new MidPointCDOEngine(yieldHandle));
        ext::shared_ptr<PricingEngine> integralCDOEngine(new IntegralCDOEngine(yieldHandle));

        const Size i = dataSet;
        correlation->setValue(hwData7[i].correlation);

// define the models to be tested
// update the tolerance values based on the models
#pragma region Models
        // std::vector<ext::shared_ptr<DefaultLossModel>> basketModels;
        // std::vector<std::string> modelNames;
        // std::vector<Real> relativeToleranceMidp, relativeTolerancePeriod, absoluteTolerance;

        std::vector<model> models;

        // safety check: if nm or nm are smaller than -1, then the model is not defined
        // return directly

        if (hwData7[i].nm < -1 || hwData7[i].nz < -1) {
            return;
        } 


        // choose distribution base on nm and nz
        if (hwData7[i].nm == -1 && hwData7[i].nz == -1) {
            ext::shared_ptr<GaussianConstantLossLM> gaussKtLossLM(
                                            new GaussianConstantLossLM(
                                                hCorrelation, 
                                                std::vector<Real>(poolSize, recovery),
                                                LatentModelIntegrationType::GaussianQuadrature, poolSize,
                                                GaussianCopulaPolicy::initTraits()
                                                )
                                           );
            
            // 1.-Inhomogeneous gaussian
            models.emplace_back(
                "Inhomogeneous gaussian",   // model name
                ext::shared_ptr<DefaultLossModel>(new IHGaussPoolLossModel(gaussKtLossLM, nBuckets, 5., -5, 15)),   // model
                1.,     // absolute tolerance
                0.04,   // relative tolerance midp
                0.04    // relative tolerance period
            );

            // 2.-homogeneous gaussian
            models.emplace_back(
                "Homogeneous gaussian",   // model name
                ext::shared_ptr<DefaultLossModel>(new HomogGaussPoolLossModel(gaussKtLossLM, nBuckets, 5., -5, 15)),   // model
                1.,     // absolute tolerance
                0.04,   // relative tolerance midp
                0.04    // relative tolerance period
            );
            
            // 3.-random default gaussian
            models.emplace_back(
                "Random default gaussian",   // model name
                ext::shared_ptr<DefaultLossModel>(new RandomDefaultLM<GaussianCopulaPolicy>(gaussKtLossLM, numSims)),   // model
                1.,     // absolute tolerance
                0.07,   // relative tolerance midp
                0.07    // relative tolerance period
            );

        } 
        else // either one is -1 and the other one is > 0, or bot are > 0
        {
            TCopulaPolicy::initTraits initTG;

            // freedom degrees. if nm is -1, then use 45 as the freedom degree which is close to gaussian
            initTG.tOrders.push_back(hwData7[i].nm > 0 ? hwData7[i].nm : 45); 
            initTG.tOrders.push_back(hwData7[i].nz > 0 ? hwData7[i].nz : 45); 

            ext::shared_ptr<TConstantLossLM> TKtLossLM(new TConstantLossLM(
                hCorrelation, std::vector<Real>(poolSize, recovery),
                LatentModelIntegrationType::GaussianQuadrature, poolSize, initTG));

            // 1.-inhomogeneous            
            models.emplace_back(
                "Inhomogeneous " +  string(hwData7[i].nm > 0 ? "student" : "gaussian") +  string("-") +  string(hwData7[i].nz > 0 ? "student" : "gaussian"),   // model name
                ext::shared_ptr<DefaultLossModel>(new IHStudentPoolLossModel(TKtLossLM, nBuckets, 5., -5., 15)),   // model
                1.,     // absolute tolerance
                0.04,   // relative tolerance midp
                0.04    // relative tolerance period
            );


            // 2.-homogeneous
            models.emplace_back(
                "Homogeneous " +  string(hwData7[i].nm > 0 ? "student" : "gaussian") +  string("-") +  string(hwData7[i].nz > 0 ? "student" : "gaussian"),   // model name
                ext::shared_ptr<DefaultLossModel>(new HomogTPoolLossModel(TKtLossLM, nBuckets, 5., -5., 15)),   // model
                1.,     // absolute tolerance
                0.04,   // relative tolerance midp
                0.04    // relative tolerance period
            );

            // 3.-random default
            models.emplace_back(
                "Random default " +  string(hwData7[i].nm > 0 ? "student" : "gaussian") +  string("-") +  string(hwData7[i].nz > 0 ? "student" : "gaussian"),   // model name
                ext::shared_ptr<DefaultLossModel>(new RandomDefaultLM<TCopulaPolicy>(TKtLossLM, numSims)),   // model
                1.,     // absolute tolerance
                0.07,   // relative tolerance midp
                0.07    // relative tolerance period
            );
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

            for (Size im = 0; im < models.size(); im++) {

                basketPtr->setLossModel(models[im].basketModel);

                cdoe.setPricingEngine(midPCDOEngine);
                check(
                        i, // data set
                        j, // tranche
                        models[im].modelName + std::string(" with midp integration on ") + trancheId.str(), // description
                        cdoe.fairPremium() * 1e4, // found 
                        hwData7[i].trancheSpread[j], // expected
                        models[im].absoluteTolerance, // absolute tolerance
                        models[im].relativeToleranceMidp // relative tolerance
                    );

                cdoe.setPricingEngine(integralCDOEngine);
                check(
                        i, // data set
                        j, // tranche
                        models[im].modelName + std::string(" with step integration on ") + trancheId.str(), // description
                        cdoe.fairPremium() * 1e4, // found
                        hwData7[i].trancheSpread[j], // expected
                        models[im].absoluteTolerance, // absolute tolerance
                        models[im].relativeTolerancePeriod // relative tolerance
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
