/*
// Copyright (c) 2015 Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#include <Hoa.hpp>
#include "catch.hpp"
#include <cfloat>
#include <cmath>
#include <iostream>

using namespace hoa;

TEST_CASE("Harmonics 2D", "[Harmonics] [2D]")
{
    typedef Harmonic<Hoa2d, float> Harmonic;
    Harmonic Harmonic0(0);
    Harmonic Harmonic1(1);
    Harmonic Harmonic2(2);
    Harmonic Harmonic3(3);
    Harmonic Harmonic4(4);
    Harmonic Harmonic5(5);
    Harmonic Harmonic6(6);
    Harmonic Harmonic13(13);
    Harmonic Harmonic14(14);
    Harmonic Harmonic21(21);
    Harmonic Harmonic22(22);

    SECTION("Static Degree")
    {
        CHECK(Harmonic::getDegree(0) == 0);
        CHECK(Harmonic::getDegree(1) == 1);
        CHECK(Harmonic::getDegree(2) == 1);
        CHECK(Harmonic::getDegree(3) == 2);
        CHECK(Harmonic::getDegree(4) == 2);
        CHECK(Harmonic::getDegree(5) == 3);
        CHECK(Harmonic::getDegree(6) == 3);
        CHECK(Harmonic::getDegree(13) == 7);
        CHECK(Harmonic::getDegree(14) == 7);
        CHECK(Harmonic::getDegree(21) == 11);
        CHECK(Harmonic::getDegree(22) == 11);
    }

    SECTION("Static Order")
    {
        CHECK(Harmonic::getOrder(0) == 0);
        CHECK(Harmonic::getOrder(1) == -1);
        CHECK(Harmonic::getOrder(2) == 1);
        CHECK(Harmonic::getOrder(3) == -2);
        CHECK(Harmonic::getOrder(4) == 2);
        CHECK(Harmonic::getOrder(5) == -3);
        CHECK(Harmonic::getOrder(6) == 3);
        CHECK(Harmonic::getOrder(13) == -7);
        CHECK(Harmonic::getOrder(14) == 7);
        CHECK(Harmonic::getOrder(21) == -11);
        CHECK(Harmonic::getOrder(22) == 11);
    }

    SECTION("Static Index")
    {
        CHECK(Harmonic::getIndex(0, 0) == 0);
        CHECK(Harmonic::getIndex(1, -1) == 1);
        CHECK(Harmonic::getIndex(1, 1) == 2);
        CHECK(Harmonic::getIndex(2, -2) == 3);
        CHECK(Harmonic::getIndex(2, 2) == 4);
        CHECK(Harmonic::getIndex(3, -3) == 5);
        CHECK(Harmonic::getIndex(3, 3) == 6);
        CHECK(Harmonic::getIndex(7, -7) == 13);
        CHECK(Harmonic::getIndex(7, 7) == 14);
        CHECK(Harmonic::getIndex(11, -11) == 21);
        CHECK(Harmonic::getIndex(11, 11) == 22);
    }

    SECTION("Static Number of Harmonics")
    {
        CHECK(Harmonic::getNumberOfHarmonics(0) == 1);
        CHECK(Harmonic::getNumberOfHarmonics(1) == 3);
        CHECK(Harmonic::getNumberOfHarmonics(2) == 5);
        CHECK(Harmonic::getNumberOfHarmonics(3) == 7);
        CHECK(Harmonic::getNumberOfHarmonics(7) == 15);
        CHECK(Harmonic::getNumberOfHarmonics(11) == 23);
    }
    
    SECTION("Static Number of Harmonics In Degree")
    {
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(0) == 1);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(1) == 2);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(2) == 2);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(3) == 2);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(7) == 2);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(11) == 2);
    }

    SECTION("Static Normalization")
    {
        CHECK(Harmonic::getSemiNormalization(0, 0) == 1.);
        CHECK(Harmonic::getSemiNormalization(1, -1) == 1.);
        CHECK(Harmonic::getSemiNormalization(1, 1) == 1.);
        CHECK(Harmonic::getSemiNormalization(2, -2) == 1.);
        CHECK(Harmonic::getSemiNormalization(2, 2) == 1.);
        CHECK(Harmonic::getSemiNormalization(3, -3) == 1.);
        CHECK(Harmonic::getSemiNormalization(3, 3) == 1.);
        CHECK(Harmonic::getSemiNormalization(7, -7) == 1.);
        CHECK(Harmonic::getSemiNormalization(7, 7) == 1.);
        CHECK(Harmonic::getSemiNormalization(11, -11) == 1.);
        CHECK(Harmonic::getSemiNormalization(11, 11) == 1.);
    }

    SECTION("Static Semi Normalization")
    {
        CHECK(Harmonic::getSemiNormalization(0, 0) == 1.);
        CHECK(Harmonic::getSemiNormalization(1, -1) == 1.);
        CHECK(Harmonic::getSemiNormalization(1, 1) == 1.);
        CHECK(Harmonic::getSemiNormalization(2, -2) == 1.);
        CHECK(Harmonic::getSemiNormalization(2, 2) == 1.);
        CHECK(Harmonic::getSemiNormalization(3, -3) == 1.);
        CHECK(Harmonic::getSemiNormalization(3, 3) == 1.);
        CHECK(Harmonic::getSemiNormalization(7, -7) == 1.);
        CHECK(Harmonic::getSemiNormalization(7, 7) == 1.);
        CHECK(Harmonic::getSemiNormalization(11, -11) == 1.);
        CHECK(Harmonic::getSemiNormalization(11, 11) == 1.);
    }







    SECTION("Local Degree")
    {
        CHECK(Harmonic0.getDegree() == 0);
        CHECK(Harmonic1.getDegree() == 1);
        CHECK(Harmonic2.getDegree() == 1);
        CHECK(Harmonic3.getDegree() == 2);
        CHECK(Harmonic4.getDegree() == 2);
        CHECK(Harmonic5.getDegree() == 3);
        CHECK(Harmonic6.getDegree() == 3);
        CHECK(Harmonic13.getDegree() == 7);
        CHECK(Harmonic14.getDegree() == 7);
        CHECK(Harmonic21.getDegree() == 11);
        CHECK(Harmonic22.getDegree() == 11);
    }

    SECTION("Local Order")
    {
        CHECK(Harmonic0.getOrder() == 0);
        CHECK(Harmonic1.getOrder() == -1);
        CHECK(Harmonic2.getOrder() == 1);
        CHECK(Harmonic3.getOrder() == -2);
        CHECK(Harmonic4.getOrder() == 2);
        CHECK(Harmonic5.getOrder() == -3);
        CHECK(Harmonic6.getOrder() == 3);
        CHECK(Harmonic13.getOrder() == -7);
        CHECK(Harmonic14.getOrder() == 7);
        CHECK(Harmonic21.getOrder() == -11);
        CHECK(Harmonic22.getOrder() == 11);
    }

    SECTION("Local Index")
    {
        CHECK(Harmonic0.getIndex() == 0);
        CHECK(Harmonic1.getIndex() == 1);
        CHECK(Harmonic2.getIndex() == 2);
        CHECK(Harmonic3.getIndex() == 3);
        CHECK(Harmonic4.getIndex() == 4);
        CHECK(Harmonic5.getIndex() == 5);
        CHECK(Harmonic6.getIndex() == 6);
        CHECK(Harmonic13.getIndex() == 13);
        CHECK(Harmonic14.getIndex() == 14);
        CHECK(Harmonic21.getIndex() == 21);
        CHECK(Harmonic22.getIndex() == 22);
    }

    SECTION("Local Normalization")
    {
        CHECK(Harmonic0.getNormalization() == 1.);
        CHECK(Harmonic1.getNormalization() == 1.);
        CHECK(Harmonic2.getNormalization() == 1.);
        CHECK(Harmonic3.getNormalization() == 1.);
        CHECK(Harmonic4.getNormalization() == 1.);
        CHECK(Harmonic5.getNormalization() == 1.);
        CHECK(Harmonic6.getNormalization() == 1.);
        CHECK(Harmonic13.getNormalization() == 1.);
        CHECK(Harmonic14.getNormalization() == 1.);
        CHECK(Harmonic21.getNormalization() == 1.);
        CHECK(Harmonic22.getNormalization() == 1.);
    }

    SECTION("Local Semi Normalization")
    {
        CHECK(Harmonic0.getNormalization() == 1.);
        CHECK(Harmonic1.getNormalization() == 1.);
        CHECK(Harmonic2.getNormalization() == 1.);
        CHECK(Harmonic3.getNormalization() == 1.);
        CHECK(Harmonic4.getNormalization() == 1.);
        CHECK(Harmonic5.getNormalization() == 1.);
        CHECK(Harmonic6.getNormalization() == 1.);
        CHECK(Harmonic13.getNormalization() == 1.);
        CHECK(Harmonic14.getNormalization() == 1.);
        CHECK(Harmonic21.getNormalization() == 1.);
        CHECK(Harmonic22.getNormalization() == 1.);
    }
}



TEST_CASE("Harmonics 3D", "[Harmonics] [3D]")
{
    const float epsilon = FLT_EPSILON;
    typedef Harmonic<Hoa3d, float> Harmonic;
    Harmonic Harmonic0(0);

    Harmonic Harmonic1(1);
    Harmonic Harmonic2(2);
    Harmonic Harmonic3(3);

    Harmonic Harmonic9(9);
    Harmonic Harmonic10(10);
    Harmonic Harmonic11(11);
    Harmonic Harmonic12(12);
    Harmonic Harmonic13(13);
    Harmonic Harmonic14(14);
    Harmonic Harmonic15(15);

    Harmonic Harmonic49(49);
    Harmonic Harmonic50(50);
    Harmonic Harmonic51(51);
    Harmonic Harmonic52(52);
    Harmonic Harmonic53(53);
    Harmonic Harmonic54(54);
    Harmonic Harmonic55(55);
    Harmonic Harmonic56(56);
    Harmonic Harmonic57(57);
    Harmonic Harmonic58(58);
    Harmonic Harmonic59(59);
    Harmonic Harmonic60(60);
    Harmonic Harmonic61(61);
    Harmonic Harmonic62(62);
    Harmonic Harmonic63(63);

    SECTION("Static Degree")
    {
        CHECK(Harmonic::getDegree(0) == 0);

        CHECK(Harmonic::getDegree(1) == 1);
        CHECK(Harmonic::getDegree(2) == 1);
        CHECK(Harmonic::getDegree(3) == 1);

        CHECK(Harmonic::getDegree(9) == 3);
        CHECK(Harmonic::getDegree(10) == 3);
        CHECK(Harmonic::getDegree(11) == 3);
        CHECK(Harmonic::getDegree(12) == 3);
        CHECK(Harmonic::getDegree(13) == 3);
        CHECK(Harmonic::getDegree(14) == 3);
        CHECK(Harmonic::getDegree(15) == 3);

        CHECK(Harmonic::getDegree(49) == 7);
        CHECK(Harmonic::getDegree(50) == 7);
        CHECK(Harmonic::getDegree(51) == 7);
        CHECK(Harmonic::getDegree(52) == 7);
        CHECK(Harmonic::getDegree(53) == 7);
        CHECK(Harmonic::getDegree(54) == 7);
        CHECK(Harmonic::getDegree(55) == 7);
        CHECK(Harmonic::getDegree(56) == 7);
        CHECK(Harmonic::getDegree(57) == 7);
        CHECK(Harmonic::getDegree(58) == 7);
        CHECK(Harmonic::getDegree(59) == 7);
        CHECK(Harmonic::getDegree(60) == 7);
        CHECK(Harmonic::getDegree(61) == 7);
        CHECK(Harmonic::getDegree(62) == 7);
        CHECK(Harmonic::getDegree(63) == 7);
    }

    SECTION("Static Order")
    {
        CHECK(Harmonic::getOrder(0) == 0);

        CHECK(Harmonic::getOrder(1) == -1);
        CHECK(Harmonic::getOrder(2) == 0);
        CHECK(Harmonic::getOrder(3) == 1);

        CHECK(Harmonic::getOrder(9) == -3);
        CHECK(Harmonic::getOrder(10) == -2);
        CHECK(Harmonic::getOrder(11) == -1);
        CHECK(Harmonic::getOrder(12) == 0);
        CHECK(Harmonic::getOrder(13) == 1);
        CHECK(Harmonic::getOrder(14) == 2);
        CHECK(Harmonic::getOrder(15) == 3);

        CHECK(Harmonic::getOrder(49) == -7);
        CHECK(Harmonic::getOrder(50) == -6);
        CHECK(Harmonic::getOrder(51) == -5);
        CHECK(Harmonic::getOrder(52) == -4);
        CHECK(Harmonic::getOrder(53) == -3);
        CHECK(Harmonic::getOrder(54) == -2);
        CHECK(Harmonic::getOrder(55) == -1);
        CHECK(Harmonic::getOrder(56) == 0);
        CHECK(Harmonic::getOrder(57) == 1);
        CHECK(Harmonic::getOrder(58) == 2);
        CHECK(Harmonic::getOrder(59) == 3);
        CHECK(Harmonic::getOrder(60) == 4);
        CHECK(Harmonic::getOrder(61) == 5);
        CHECK(Harmonic::getOrder(62) == 6);
        CHECK(Harmonic::getOrder(63) == 7);
    }

    SECTION("Static Index")
    {
        CHECK(Harmonic::getIndex(0, 0) == 0);

        CHECK(Harmonic::getIndex(1, -1) == 1);
        CHECK(Harmonic::getIndex(1, 0) == 2);
        CHECK(Harmonic::getIndex(1, 1) == 3);

        CHECK(Harmonic::getIndex(3, -3) == 9);
        CHECK(Harmonic::getIndex(3, -2) == 10);
        CHECK(Harmonic::getIndex(3, -1) == 11);
        CHECK(Harmonic::getIndex(3, 0) == 12);
        CHECK(Harmonic::getIndex(3, 1) == 13);
        CHECK(Harmonic::getIndex(3, 2) == 14);
        CHECK(Harmonic::getIndex(3, 3) == 15);

        CHECK(Harmonic::getIndex(7, -7) == 49);
        CHECK(Harmonic::getIndex(7, -6) == 50);
        CHECK(Harmonic::getIndex(7, -5) == 51);
        CHECK(Harmonic::getIndex(7, -4) == 52);
        CHECK(Harmonic::getIndex(7, -3) == 53);
        CHECK(Harmonic::getIndex(7, -2) == 54);
        CHECK(Harmonic::getIndex(7, -1) == 55);
        CHECK(Harmonic::getIndex(7, 0) == 56);
        CHECK(Harmonic::getIndex(7, 1) == 57);
        CHECK(Harmonic::getIndex(7, 2) == 58);
        CHECK(Harmonic::getIndex(7, 3) == 59);
        CHECK(Harmonic::getIndex(7, 4) == 60);
        CHECK(Harmonic::getIndex(7, 5) == 61);
        CHECK(Harmonic::getIndex(7, 6) == 62);
        CHECK(Harmonic::getIndex(7, 7) == 63);
    }


    SECTION("Static Number of Harmonics")
    {
        CHECK(Harmonic::getNumberOfHarmonics(0) == 1);
        CHECK(Harmonic::getNumberOfHarmonics(1) == 4);
        CHECK(Harmonic::getNumberOfHarmonics(2) == 9);
        CHECK(Harmonic::getNumberOfHarmonics(3) == 16);
        CHECK(Harmonic::getNumberOfHarmonics(7) == 64);
        CHECK(Harmonic::getNumberOfHarmonics(11) == 144);
    }
    
    SECTION("Static Number of Harmonics In Degree")
    {
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(0) == 1);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(1) == 3);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(2) == 5);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(3) == 7);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(7) == 15);
        CHECK(Harmonic::getNumberOfHarmonicsInDegree(11) == 23);
    }

    
    SECTION("Static Normalization")
    {
        CHECK(Harmonic::getNormalization(0, 0) == 1.);

        CHECK(std::abs(Harmonic::getNormalization(1, -1) - 1.7320507764816284) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(1, 0) -  1.7320507764816284) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(1, 1) -  1.7320507764816284) < epsilon);
        
        CHECK(std::abs(Harmonic::getNormalization(3, -3) - 0.13944333791732788) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(3, -2) - 0.34156504273414612) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(3, -1) - 1.0801234245300293) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(3, 0) - 2.6457512378692627) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(3, 1) - 1.0801234245300293) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(3, 2) - 0.34156504273414612) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(3, 3) - 0.13944333791732788) < epsilon);
        
    
        CHECK(std::abs(Harmonic::getNormalization(7, -7) - 0.000018550535969552584) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, -6) - 0.000069409754360094666) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, -5) - 0.00035392167046666145) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, -4) - 0.0021235300227999687) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, -3) - 0.014085904695093632) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, -2) - 0.099602386355400085) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, -1) - 0.73192507028579712) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, 0) - 3.872983455657959) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, 1) - 0.73192507028579712) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, 2) - 0.099602386355400085) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, 3) - 0.014085904695093632) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, 4) - 0.0021235300227999687) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, 5) - 0.00035392167046666145) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, 6) - 0.000069409754360094666) < epsilon);
        CHECK(std::abs(Harmonic::getNormalization(7, 7) - 0.000018550535969552584) < epsilon);
    }
    
    SECTION("Static Semi Normalization")
    {
        CHECK(Harmonic::getSemiNormalization(0, 0) == 1.);
        
        CHECK(std::abs(Harmonic::getSemiNormalization(1, -1) - 1.) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(1, 0) - 1) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(1, 1) - 1.) < epsilon);
        
        CHECK(std::abs(Harmonic::getSemiNormalization(3, -3) - 0.05270462765) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(3, -2) - 0.1290994448) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(3, -1) - 0.4082482904) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(3, 0) - 1) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(3, 1) - 0.4082482904) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(3, 2) - 0.1290994448) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(3, 3) - 0.05270462765) < epsilon);
        
        CHECK(std::abs(Harmonic::getSemiNormalization(7, -7) - 0.000004789727674) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, -6) - 0.00001792151992) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, -5) - 0.00009138217984) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, -4) - 0.0005482930789) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, -3) - 0.003636964836) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, -2) - 0.02571722499) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, -1) - 0.1889822365) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, 0) - 1) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, 1) - 0.1889822365) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, 2) - 0.02571722499) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, 3) - 0.003636964836) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, 4) - 0.0005482930789) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, 5) - 0.00009138217984) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, 6) - 0.00001792151992) < epsilon);
        CHECK(std::abs(Harmonic::getSemiNormalization(7, 7) - 0.000004789727674) < epsilon);
    }





    SECTION("Local Degree")
    {
        CHECK(Harmonic0.getDegree() == 0);
        
        CHECK(Harmonic1.getDegree() == 1);
        CHECK(Harmonic2.getDegree() == 1);
        CHECK(Harmonic3.getDegree() == 1);
        
        CHECK(Harmonic9.getDegree() == 3);
        CHECK(Harmonic10.getDegree() == 3);
        CHECK(Harmonic11.getDegree() == 3);
        CHECK(Harmonic12.getDegree() == 3);
        CHECK(Harmonic13.getDegree() == 3);
        CHECK(Harmonic14.getDegree() == 3);
        CHECK(Harmonic15.getDegree() == 3);
        
        CHECK(Harmonic49.getDegree() == 7);
        CHECK(Harmonic50.getDegree() == 7);
        CHECK(Harmonic51.getDegree() == 7);
        CHECK(Harmonic52.getDegree() == 7);
        CHECK(Harmonic53.getDegree() == 7);
        CHECK(Harmonic54.getDegree() == 7);
        CHECK(Harmonic55.getDegree() == 7);
        CHECK(Harmonic56.getDegree() == 7);
        CHECK(Harmonic57.getDegree() == 7);
        CHECK(Harmonic58.getDegree() == 7);
        CHECK(Harmonic59.getDegree() == 7);
        CHECK(Harmonic60.getDegree() == 7);
        CHECK(Harmonic61.getDegree() == 7);
        CHECK(Harmonic62.getDegree() == 7);
        CHECK(Harmonic63.getDegree() == 7);
    }

    SECTION("Local Order")
    {
        CHECK(Harmonic0.getOrder()== 0);
        
        CHECK(Harmonic1.getOrder()== -1);
        CHECK(Harmonic2.getOrder()== 0);
        CHECK(Harmonic3.getOrder()== 1);
        
        CHECK(Harmonic9.getOrder()== -3);
        CHECK(Harmonic10.getOrder()== -2);
        CHECK(Harmonic11.getOrder()== -1);
        CHECK(Harmonic12.getOrder()== 0);
        CHECK(Harmonic13.getOrder()== 1);
        CHECK(Harmonic14.getOrder()== 2);
        CHECK(Harmonic15.getOrder()== 3);
        
        CHECK(Harmonic49.getOrder()== -7);
        CHECK(Harmonic50.getOrder()== -6);
        CHECK(Harmonic51.getOrder()== -5);
        CHECK(Harmonic52.getOrder()== -4);
        CHECK(Harmonic53.getOrder()== -3);
        CHECK(Harmonic54.getOrder()== -2);
        CHECK(Harmonic55.getOrder()== -1);
        CHECK(Harmonic56.getOrder()== 0);
        CHECK(Harmonic57.getOrder()== 1);
        CHECK(Harmonic58.getOrder()== 2);
        CHECK(Harmonic59.getOrder()== 3);
        CHECK(Harmonic60.getOrder()== 4);
        CHECK(Harmonic61.getOrder()== 5);
        CHECK(Harmonic62.getOrder()== 6);
        CHECK(Harmonic63.getOrder()== 7);
    }
    
    SECTION("Local Index")
    {
        CHECK(Harmonic0.getIndex()== 0);
        
        CHECK(Harmonic1.getIndex()== 1);
        CHECK(Harmonic2.getIndex()== 2);
        CHECK(Harmonic3.getIndex()== 3);
        
        CHECK(Harmonic9.getIndex()== 9);
        CHECK(Harmonic10.getIndex()== 10);
        CHECK(Harmonic11.getIndex()== 11);
        CHECK(Harmonic12.getIndex()== 12);
        CHECK(Harmonic13.getIndex()== 13);
        CHECK(Harmonic14.getIndex()== 14);
        CHECK(Harmonic15.getIndex()== 15);
        
        CHECK(Harmonic49.getIndex()== 49);
        CHECK(Harmonic50.getIndex()== 50);
        CHECK(Harmonic51.getIndex()== 51);
        CHECK(Harmonic52.getIndex()== 52);
        CHECK(Harmonic53.getIndex()== 53);
        CHECK(Harmonic54.getIndex()== 54);
        CHECK(Harmonic55.getIndex()== 55);
        CHECK(Harmonic56.getIndex()== 56);
        CHECK(Harmonic57.getIndex()== 57);
        CHECK(Harmonic58.getIndex()== 58);
        CHECK(Harmonic59.getIndex()== 59);
        CHECK(Harmonic60.getIndex()== 60);
        CHECK(Harmonic61.getIndex()== 61);
        CHECK(Harmonic62.getIndex()== 62);
        CHECK(Harmonic63.getIndex()== 63);
    }

    SECTION("Local Normalization")
    {
        CHECK(Harmonic0.getNormalization() == 1.);
        
        CHECK(std::abs(Harmonic1.getNormalization() - 1.7320507764816284) < epsilon);
        CHECK(std::abs(Harmonic2.getNormalization() - 1.7320507764816284) < epsilon);
        CHECK(std::abs(Harmonic3.getNormalization() - 1.7320507764816284) < epsilon);
        
        CHECK(std::abs(Harmonic9.getNormalization() - 0.13944333791732788) < epsilon);
        CHECK(std::abs(Harmonic10.getNormalization() - 0.34156504273414612) < epsilon);
        CHECK(std::abs(Harmonic11.getNormalization() - 1.0801234245300293) < epsilon);
        CHECK(std::abs(Harmonic12.getNormalization() - 2.6457512378692627) < epsilon);
        CHECK(std::abs(Harmonic13.getNormalization() - 1.0801234245300293) < epsilon);
        CHECK(std::abs(Harmonic14.getNormalization() - 0.34156504273414612) < epsilon);
        CHECK(std::abs(Harmonic15.getNormalization() - 0.13944333791732788) < epsilon);
        
        CHECK(std::abs(Harmonic49.getNormalization() - 0.000018550535969552584) < epsilon);
        CHECK(std::abs(Harmonic50.getNormalization() - 0.000069409754360094666) < epsilon);
        CHECK(std::abs(Harmonic51.getNormalization() - 0.00035392167046666145) < epsilon);
        CHECK(std::abs(Harmonic52.getNormalization() - 0.0021235300227999687) < epsilon);
        CHECK(std::abs(Harmonic53.getNormalization() - 0.014085904695093632) < epsilon);
        CHECK(std::abs(Harmonic54.getNormalization() - 0.099602386355400085) < epsilon);
        CHECK(std::abs(Harmonic55.getNormalization() - 0.73192507028579712) < epsilon);
        CHECK(std::abs(Harmonic56.getNormalization() - 1) < 3.872983455657959);
        CHECK(std::abs(Harmonic57.getNormalization() - 0.73192507028579712) < epsilon);
        CHECK(std::abs(Harmonic58.getNormalization() - 0.099602386355400085) < epsilon);
        CHECK(std::abs(Harmonic59.getNormalization() - 0.014085904695093632) < epsilon);
        CHECK(std::abs(Harmonic60.getNormalization() - 0.0021235300227999687) < epsilon);
        CHECK(std::abs(Harmonic61.getNormalization() - 0.00035392167046666145) < epsilon);
        CHECK(std::abs(Harmonic62.getNormalization() - 0.000069409754360094666) < epsilon);
        CHECK(std::abs(Harmonic63.getNormalization() - 0.000018550535969552584) < epsilon);
    }
    
    SECTION("Local Semi Normalization")
    {
        CHECK(Harmonic0.getSemiNormalization() == 1.);
        
        CHECK(std::abs(Harmonic1.getSemiNormalization() - 1.) < epsilon);
        CHECK(std::abs(Harmonic2.getSemiNormalization() - 1) < epsilon);
        CHECK(std::abs(Harmonic3.getSemiNormalization() - 1.) < epsilon);
        
        CHECK(std::abs(Harmonic9.getSemiNormalization() - 0.05270462765) < epsilon);
        CHECK(std::abs(Harmonic10.getSemiNormalization() - 0.1290994448) < epsilon);
        CHECK(std::abs(Harmonic11.getSemiNormalization() - 0.4082482904) < epsilon);
        CHECK(std::abs(Harmonic12.getSemiNormalization() - 1) < epsilon);
        CHECK(std::abs(Harmonic13.getSemiNormalization() - 0.4082482904) < epsilon);
        CHECK(std::abs(Harmonic14.getSemiNormalization() - 0.1290994448) < epsilon);
        CHECK(std::abs(Harmonic15.getSemiNormalization() - 0.05270462765) < epsilon);
        
        CHECK(std::abs(Harmonic49.getSemiNormalization() - 0.000004789727674) < epsilon);
        CHECK(std::abs(Harmonic50.getSemiNormalization() - 0.00001792151992) < epsilon);
        CHECK(std::abs(Harmonic51.getSemiNormalization() - 0.00009138217984) < epsilon);
        CHECK(std::abs(Harmonic52.getSemiNormalization() - 0.0005482930789) < epsilon);
        CHECK(std::abs(Harmonic53.getSemiNormalization() - 0.003636964836) < epsilon);
        CHECK(std::abs(Harmonic54.getSemiNormalization() - 0.02571722499) < epsilon);
        CHECK(std::abs(Harmonic55.getSemiNormalization() - 0.1889822365) < epsilon);
        CHECK(std::abs(Harmonic56.getSemiNormalization() - 1) < epsilon);
        CHECK(std::abs(Harmonic57.getSemiNormalization() - 0.1889822365) < epsilon);
        CHECK(std::abs(Harmonic58.getSemiNormalization() - 0.02571722499) < epsilon);
        CHECK(std::abs(Harmonic59.getSemiNormalization() - 0.003636964836) < epsilon);
        CHECK(std::abs(Harmonic60.getSemiNormalization() - 0.0005482930789) < epsilon);
        CHECK(std::abs(Harmonic61.getSemiNormalization() - 0.00009138217984) < epsilon);
        CHECK(std::abs(Harmonic62.getSemiNormalization() - 0.00001792151992) < epsilon);
        CHECK(std::abs(Harmonic63.getSemiNormalization() - 0.000004789727674) < epsilon);
    }
}
