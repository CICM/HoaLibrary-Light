/*
// Copyright (c) 2012-2017 CICM - Universite Paris 8 - Labex Arts H2H.
// Authors :
// 2012: Pierre Guillot, Eliott Paris & Julien Colafrancesco.
// 2012-2015: Pierre Guillot & Eliott Paris.
// 2015: Pierre Guillot & Eliott Paris & Thomas Le Meur (Light version)
// 2016-2017: Pierre Guillot.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#pragma once

#include "Hoa_Processor.hpp"

namespace hoa
{
    //! The echanger class renumber and normalize the harmonics channels.
    /** The echanger should be used to renumber and normalize the harmonics channels. The library uses the Ambisonics Channels Numbering (ACN), this class allows to convert channels arrengements from Furse-Malham (B-format) or Single Index (SID) to  Ambisonics Channels Numbering (ACN) and conversely. Furse-Malham and SID never reach up to 3rd order so the maximum order of decomposition should be 3. The library uses the semi-normalization (SN2D and SN3D), this class allows to normalize the channels to the full normalization (N2D and N3D) or to MaxN (B-format) and conversely.
     */
    template <Dimension D, typename T> class Exchanger : public ProcessorHarmonics<D, T>
    {
    public:

        //! The numbering conversion.
        /** The enum defines the numbering conversion.
         */
        enum Numbering
        {
            ACN             = 0, /*!<  The numbering is considered as ACN. */
            fromFurseMalham = 1, /*!<  From Furse-Malham (B-format) to ACN. */
            fromSID         = 2, /*!<  From SID to ACN. */
            toFurseMalham   = 3, /*!<  To Furse-Malham (B-format) from ACN. */
            toSID           = 4  /*!<  To SID from ACN. */
        };

        //! The numbering conversion.
        /** The enum defines the numbering conversion.
         */
        enum Normalization
        {
            SN2D            = 0, /*!<  The normalization is considered as semi-normalization. */
            SN3D            = 0, /*!<  The normalization is considered as semi-normalization. */
            fromN3D         = 0, /*!<  From N3D to SN3D. */
            fromMaxN        = 1, /*!<  From MaxN (B-format) to SN2D/SN3D. */
            toN3D           = 2, /*!<  To N3D from SN3D. */
            toMaxN          = 3  /*!<  To MaxN (B-format) from SN2D/SN3D. */
        };

        //! The exchanger constructor.
        /**	The exchanger constructor allocates and initialize the member values to renumber and normalize the harmonics channels. The order must be at least 1 and should be 3 at maximum.
         @param     order	The order.
         */
        Exchanger(const size_t order) noexcept;

        //! The exchanger destructor.
        /**	The exchanger destructor free the memory.
         */
		virtual ~Exchanger() noexcept = 0;

        //! This method performs the numbering and the normalization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        virtual void process(T const* inputs, T* outputs) noexcept;

        //! Sets the numbering conversion.
        /**	This method sets the numbering conversion.
         @param mode The numbering convertion.
         */
        virtual void setNumbering(const Numbering mode) noexcept;

        //! Gets the numbering conversion.
        /**	This method gets the numbering conversion.
         @return The numbering convertion.
         */
        virtual Numbering getNumbering() const noexcept;
        
        //! Sets the normalization conversion.
        /**	This method sets the normalization conversion.
         @param mode The normalization convertion.
         */
        virtual void setNormalization(const Normalization mode) noexcept;
        
        //! Gets the normalization conversion.
        /**	This method gets the normalization conversion.
         @return The normalization convertion.
         */
        virtual Normalization getNormalization() const noexcept;

        //! Retrieves the harmonic order of an input depending on the current numbering configuration.
        /**
         @param     index	The index of an harmonic.
         @return    The harmonic order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        virtual long getInputHarmonicOrder(const size_t index) const noexcept;

        //! Retrieves the harmonic order of an output depending on the current numbering configuration.
        /**
         @param     index	The index of an harmonic.
         @return    The harmonic order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        virtual long getOutputHarmonicOrder(const size_t index) const noexcept;

        //! Retrieves the name of an harmonic depending on the current numbering configuration.
        /** This methods returns the name that contains the degree and the order of the harmonic for ACN and SID or the letter code of the harmonic for Furse-Malham.
         @param     index	The index of an harmonic.
         @param     isInput	Pass true to retrieve the input harmonic name, false for an output.
         @return    The name of the harmonic.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        virtual std::string getHarmonicName(const size_t index, const bool isInput) const noexcept;
    };

//! @cond

    template <typename T> class Exchanger<Hoa2d, T> : public ProcessorHarmonics<Hoa2d, T>
    {
    public:
        //! The numbering conversion.
        /** The enum defines the numbering conversion.
         */
        enum Numbering
        {
            ACN             = 0, /*!<  The numbering is considered as ACN. */
            fromFurseMalham = 1, /*!<  From Furse-Malham (B-format) to ACN. */
            fromSID         = 2, /*!<  From SID to ACN. */
            toFurseMalham   = 3, /*!<  To Furse-Malham (B-format) from ACN. */
            toSID           = 4  /*!<  To SID from ACN. */
        };

        //! The numbering conversion.
        /** The enum defines the numbering conversion.
         */
        enum Normalization
        {
            SN2D            = 0, /*!<  The normalization is considered as semi-normalization. */
            fromMaxN        = 1, /*!<  From MaxN (B-format) to SN2D. */
            toMaxN          = 3  /*!<  To MaxN (B-format) from SN2D. */
        };

    private:
        Numbering       m_numbering;
        Normalization   m_normalization;
    public:

        //! The exchanger constructor.
        /**	The exchanger constructor allocates and initialize the member values to renumber and normalize the harmonics channels. The order must be at least 1 and should be 3 at maximum.
         @param     order	The order.
         */
        inline Exchanger(const size_t order) noexcept : ProcessorHarmonics<Hoa2d, T>(order),
        m_numbering(ACN),
        m_normalization(SN2D)
        {
            ;
        }

        //! The exchanger destructor.
        /**	The exchanger destructor free the memory.
         */
        inline ~Exchanger() noexcept
        {
            ;
        }

        //! Sets the numbering and the normalization conversion from B-Format.
        /**	This method the numbering and the normalization conversion from B-Format. Similar to from Furse-Malham numebring and from MaxN normalization.
         */
        inline void setFromBFormat() noexcept
        {
            m_numbering = fromFurseMalham;
            m_normalization = fromMaxN;
        }

        //! Sets the numbering and the normalization conversion to B-Format.
        /**	This method the numbering and the normalization conversion from B-Format. Similar to to Furse-Malham numebring and to MaxN normalization.
         */
        inline void setToBFormat() noexcept
        {
            m_numbering = toFurseMalham;
            m_normalization = toMaxN;
        }

        //! Sets the numbering conversion.
        /**	This method sets the numbering conversion.
         @param mode The numbering convertion.
         */
        inline void setNumbering(const Numbering mode) noexcept
        {
            m_numbering = mode;
        }

        //! Gets the numbering conversion.
        /**	This method gets the numbering conversion.
         @return The numbering convertion.
         */
        inline Numbering getNumbering() const noexcept
        {
            return m_numbering;
        }

        //! Sets the normalization conversion.
        /**	This method sets the normalization conversion.
         @param mode The normalization convertion.
         */
        inline void setNormalization(const Normalization mode) noexcept
        {
            m_normalization = mode;
        }

        //! Gets the normalization conversion.
        /**	This method gets the normalization conversion.
         @return The normalization convertion.
         */
        inline Normalization getNormalization() const noexcept
        {
            return m_normalization;
        }


        //! This method performs the numbering and the normalization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void process(T const* inputs, T* outputs) noexcept
        {
            T const* ins = inputs;
            if(m_numbering == fromFurseMalham)
            {
                numberFromFurseMalham(inputs, outputs);
                ins = outputs;
            }
            else if(m_numbering == fromSID)
            {
                numberFromSID(inputs, outputs);
                ins = outputs;
            }
            switch(m_normalization)
            {
                case fromMaxN:
                    outputs[0] = ins[0] * T(sqrt(2.));
                    Signal<T>::copy(ProcessorHarmonics<Hoa2d, T>::getNumberOfHarmonics() - 1, ins+1, outputs+1);
                    ins = outputs;
                    break;
                case toMaxN:
                    outputs[0] = ins[0] / T(sqrt(2.));
                    Signal<T>::copy(ProcessorHarmonics<Hoa2d, T>::getNumberOfHarmonics() - 1, ins+1, outputs+1);
                    ins = outputs;
                    break;
                default:
                    break;
            }
            if(m_numbering == toFurseMalham)
            {
                numberToFurseMalham(ins, outputs);
            }
            else if(m_numbering == toSID)
            {
                numberToSID(ins, outputs);
            }
        }

        //! This method number the channels from Furse-Malham to ACN.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void numberFromFurseMalham(T const* inputs, T* outputs) noexcept
        {
            T temp = inputs[1];
            *(outputs++) = inputs[0]; // W -> 0
            *(outputs++) = inputs[2]; // Y -> 1
            *(outputs++) = temp;      // X -> 2
            if(ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder() > 1ul)
            {
                temp       = inputs[3];
                *(outputs++) = inputs[4]; // V -> 3
                *(outputs++) = temp;      // U -> 4
                if(ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder() > 2ul)
                {
                    temp       = inputs[5];
                    *(outputs++) = inputs[6]; // Q -> 5
                    *(outputs++) = temp;      // U -> 6
                }
            }
        }

        //! This method number the channels from SID to ACN.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void numberFromSID(T const* inputs, T* outputs) noexcept
        {
            T temp = inputs[1];
            *(outputs++) = inputs[0]; // 0 -> 0
            *(outputs++) = inputs[2]; // 2 -> 1
            *(outputs++) = temp;      // 1 -> 2
            for(size_t i = 2; i <= ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder(); i++)
            {
                temp         = inputs[(i-1)*2+1];
                *(outputs++) = inputs[(i-1)*2+2];
                *(outputs++) = temp;
            }
        }

        //! This method number the channels from ACN to Furse-Malham.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void numberToFurseMalham(T const* inputs, T* outputs) noexcept
        {
            T temp = inputs[1];
            *(outputs++) = inputs[0]; // 0 -> W
            *(outputs++) = inputs[2]; // 2 -> X
            *(outputs++) = temp;      // 1 -> Y
            if(ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder() > 1ul)
            {
                temp       = inputs[3];
                *(outputs++) = inputs[4]; // 4 -> U
                *(outputs++) = temp;      // 3 -> V
                if(ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder() > 2ul)
                {
                    temp       = inputs[5];
                    *(outputs++) = inputs[6]; // 6 -> U
                    *(outputs++) = temp;      // 5 -> Q
                }
            }
        }

        //! This method number the channels from ACN to SID.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void numberToSID(T const* inputs, T* outputs) noexcept
        {
            T temp = inputs[1];
            *(outputs++) = inputs[0]; // 0 -> 0
            *(outputs++) = inputs[2]; // 2 -> 1
            *(outputs++) = temp;      // 1 -> 2
            for(size_t i = 2; i <= ProcessorHarmonics<Hoa2d, T>::getDecompositionOrder(); i++)
            {
                temp         = inputs[(i-1)*2+1];
                *(outputs++) = inputs[(i-1)*2+2];
                *(outputs++) = temp;
            }
        }
        
        //! Retrieve the harmonic order of an input depending on the current numbering configuration.
        /** Retrieve the harmonic order of an input depending on the current numbering configuration.
         @param     index	The index of an harmonic.
         @return    The harmonic order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        long getInputHarmonicOrder(const size_t index) const noexcept
        {
            const Numbering numb = getNumbering();
            const bool acn = (numb == ACN) || (numb == toFurseMalham) || (numb == toSID) || (ProcessorHarmonics<Hoa2d, T>::getHarmonicDegree(index) > 3);
            return ProcessorHarmonics<Hoa2d, T>::getHarmonicOrder(index) * (acn ? 1l : -1l);
        }
        
        //! Retrieve the harmonic order of an output depending on the current numbering configuration.
        /** Retrieve the harmonic order of an output depending on the current numbering configuration.
         @param     index	The index of an harmonic.
         @return    The harmonic order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        long getOutputHarmonicOrder(const size_t index) const noexcept
        {
            const Numbering numb = getNumbering();
            const bool acn = (numb == ACN) || (numb == fromFurseMalham) || (numb == fromSID) || (ProcessorHarmonics<Hoa2d, T>::getHarmonicDegree(index) > 3);
            return ProcessorHarmonics<Hoa2d, T>::getHarmonicOrder(index) * (acn ? 1l : -1l);
        }
        
        //! Retrieve the name of an harmonic depending on the current numbering configuration.
        /** Retrieve the name of an harmonic depending on the current numbering configuration.
         This returns the name of the harmonic that contains its degree and its order for ACN and SID and harmonic name with letter code for FUMA.
         @param     index	The index of an harmonic.
         @param     isInput	Pass true to retrieve the input harmonic name, false for an output.
         @return    The name of the harmonic.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        std::string getHarmonicName(const size_t index, const bool isInput) const noexcept
        {
            const Numbering numb = getNumbering();
            const bool acn = (numb == ACN) || (isInput && numb == toFurseMalham) || (isInput && numb == toSID) || (!isInput && numb == fromFurseMalham) || (!isInput && numb == fromSID);
            const bool malham = !acn && ((isInput && numb == fromFurseMalham) || (!isInput && numb == toFurseMalham));
            
            if(acn || (malham && index > 6))
            {
                // [0, 0], [1, -1], [1, 1], [2, -2], [2, 2], [3, -3], [3, 3]...
                return ProcessorHarmonics<Hoa2d, T>::getHarmonicName(index);
            }
            else if(malham)
            {
                if(index < 7)
                {
                    static const char FurseMalham_2D_Letterings[7] =
                    {
                        'W', 'X', 'Y', 'U', 'V', 'P', 'Q'
                    };
                    
                    std::string name = "Harmonic ";
                    return name += FurseMalham_2D_Letterings[index];
                }
                
                return "?";
            }
            else
            {
                // [0, 0], [1, 1], [1, -1], [2, 2], [2, -2], [3, 3], [3, -3]...
                size_t hdegree = ProcessorHarmonics<Hoa2d, T>::getHarmonicDegree(index);
                size_t horder  = ProcessorHarmonics<Hoa2d, T>::getHarmonicOrder(index) * -1l;
                std::ostringstream ostr;
                ostr <<  "Harmonic " << hdegree << " " << horder;
                return ostr.str();
            }
        }
    };
    
    
    

    template <typename T> class Exchanger<Hoa3d, T> : public ProcessorHarmonics<Hoa3d, T>
    {
    public:
        //! The numbering conversion.
        /** The enum defines the numbering conversion.
         */
        enum Numbering
        {
            ACN             = 0, /*!<  The numbering is considered as ACN. */
            fromFurseMalham = 1, /*!<  From Furse-Malham (B-format) to ACN. */
            fromSID         = 2, /*!<  From SID to ACN. */
            toFurseMalham   = 3, /*!<  To Furse-Malham (B-format) from ACN. */
            toSID           = 4  /*!<  To SID from ACN. */
        };

        //! The numbering conversion.
        /** The enum defines the numbering conversion.
         */
        enum Normalization
        {
            SN3D            = 0, /*!<  The normalization is considered as semi-normalization. */
            fromN3D         = 1, /*!<  From N3D to SN2D. */
            fromMaxN        = 2, /*!<  From MaxN (B-format) to SN3D. */
            toN3D           = 3, /*!<  To N3D from SN3D. */
            toMaxN          = 4  /*!<  To MaxN (B-format) from SN3D. */
        };

    private:

        Numbering       m_numbering;
        Normalization   m_normalization;
        T*              m_harmonics;
    public:

        //! The exchanger constructor.
        /**	The exchanger constructor allocates and initialize the member values to renumber and normalize the harmonics channels. The order must be at least 1 and should be 3 at maximum.
         @param     order	The order.
         */
        inline Exchanger(const size_t order) noexcept : ProcessorHarmonics<Hoa3d, T>(order),
        m_numbering(ACN),
        m_normalization(SN3D)
        {
            m_harmonics = Signal<T>::alloc(order*2+1);
        }

        //! The exchanger destructor.
        /**	The exchanger destructor free the memory.
         */
        inline ~Exchanger() noexcept
        {
            Signal<T>::free(m_harmonics);
        }

        //! Sets the numbering and the normalization conversion from B-Format.
        /**	This method the numbering and the normalization conversion from B-Format. Similar to from Furse-Malham numebring and from MaxN normalization.
         */
        inline void setFromBFormat() noexcept
        {
            m_numbering = fromFurseMalham;
            m_normalization = fromMaxN;
        }

        //! Sets the numbering and the normalization conversion to B-Format.
        /**	This method the numbering and the normalization conversion from B-Format. Similar to to Furse-Malham numebring and to MaxN normalization.
         */
        inline void setToBFormat() noexcept
        {
            m_numbering = toFurseMalham;
            m_normalization = toMaxN;
        }

        //! Sets the numbering and the normalization conversion from B-Format.
        /**	This method the numbering and the normalization conversion from B-Format. Similar to from Furse-Malham numebring and from MaxN normalization.
         */
        inline void setFromDaniel() noexcept
        {
            m_numbering = fromSID;
            m_normalization = fromN3D;
        }

        //! Sets the numbering and the normalization conversion to B-Format.
        /**	This method the numbering and the normalization conversion from B-Format. Similar to to Furse-Malham numebring and to MaxN normalization.
         */
        inline void setToDaniel() noexcept
        {
            m_numbering = toSID;
            m_normalization = toN3D;
        }

        //! Sets the numbering conversion.
        /**	This method sets the numbering conversion.
         @param mode The numbering convertion.
         */
        inline void setNumbering(const Numbering mode) noexcept
        {
            m_numbering = mode;
        }

        //! Gets the numbering conversion.
        /**	This method gets the numbering conversion.
         @return The numbering convertion.
         */
        inline Numbering getNumbering() const noexcept
        {
            return m_numbering;
        }

        //! Sets the normalization conversion.
        /**	This method sets the normalization conversion.
         @param mode The normalization convertion.
         */
        inline void setNormalization(const Normalization mode) noexcept
        {
            m_normalization = mode;
        }

        //! Gets the normalization conversion.
        /**	This method gets the normalization conversion.
         @return The normalization convertion.
         */
        inline Normalization getNormalization() const noexcept
        {
            return m_normalization;
        }

        //! This method performs the numbering and the normalization.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void process(T const* inputs, T* outputs) noexcept
        {
            T const* ins = inputs;
            if(m_numbering == fromFurseMalham)
            {
                numberFromFurseMalham(inputs, outputs);
                ins = outputs;
            }
            else if(m_numbering == fromSID)
            {
                numberFromSID(inputs, outputs);
                ins = outputs;
            }
            switch(m_normalization)
            {
                case fromMaxN:
                    normalizeFromMaxN(ins, outputs);
                    ins = outputs;
                    break;
                case fromN3D:
                    normalizeFromN3D(ins, outputs);
                    ins = outputs;
                    break;
                case toMaxN:
                    normalizeToMaxN(ins, outputs);
                    ins = outputs;
                    break;
                case toN3D:
                    normalizeToN3D(ins, outputs);
                    ins = outputs;
                    break;
                default:
                    break;
            }
            if(m_numbering == toFurseMalham)
            {
                numberToFurseMalham(ins, outputs);
            }
            else if(m_numbering == toSID)
            {
                numberToSID(ins, outputs);
            }
        }

        //! This method number the channels from Furse-Malham to ACN.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void numberFromFurseMalham(T const* inputs, T* outputs) noexcept
        {
            T temp = inputs[1];
            outputs[0] = inputs[0]; // W(0) -> 0
            outputs[1] = inputs[2]; // Y(2) -> 1
            outputs[2] = inputs[3]; // Z(3) -> 2
            outputs[3] = temp;      // X(1) -> 3
            if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 1ul)
            {
                temp         = inputs[4];
                T temp2      = inputs[5];
                T temp3      = inputs[7];
                outputs[4] = inputs[8]; // V(8) -> 4
                outputs[5] = inputs[6]; // T(6) -> 5
                outputs[6] = temp;      // R(4) -> 6
                outputs[7] = temp2;     // S(5) -> 7
                outputs[8] = temp3;     // U(7) -> 8
                if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 2ul)
                {
                    temp      = inputs[9];
                    temp2     = inputs[10];
                    temp3     = inputs[12];
                    const T temp4   = inputs[14];
                    outputs[9]  = inputs[15];// Q(15) -> 9
                    outputs[10] = inputs[13];// O(13) -> 10
                    outputs[11] = inputs[11];// M(11) -> 11
                    outputs[12] = temp;      // K(9)  -> 12
                    outputs[13] = temp2;     // L(10) -> 13
                    outputs[14] = temp3;     // N(12) -> 14
                    outputs[15] = temp4;     // P(14) -> 15
                    if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 3ul)
                    {
                        Signal<T>::copy(ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics() - 16ul, inputs+16, outputs+16);
                    }
                }
            }
        }

        //! This method number the channels from ACN to Furse-Malham.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void numberToFurseMalham(T const* inputs, T* outputs) noexcept
        {
            T temp  = inputs[1];
            T temp2 = inputs[2];
            outputs[0] = inputs[0]; // 0 -> W(0)
            outputs[1] = inputs[3]; // 3 -> X(1)
            outputs[2] = temp;      // 1 -> Y(2)
            outputs[3] = temp2;     // 2 -> Z(3)
            if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 1ul)
            {
                temp       = inputs[4];
                temp2      = inputs[5];
                outputs[4] = inputs[6]; // 6 -> R(4)
                outputs[5] = inputs[7]; // 7 -> S(5)
                outputs[6] = temp2;     // 5 -> T(6)
                outputs[7] = inputs[8]; // 8 -> U(7)
                outputs[8] = temp;      // 4 -> V(8)
                if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 2ul)
                {
                    temp            = inputs[9];
                    temp2           = inputs[10];
                    outputs[9]  = inputs[12];// 12 -> K(9)
                    outputs[10] = inputs[13];// 13 -> L(10)
                    outputs[11] = inputs[11];// 11 -> M(11)
                    outputs[12] = inputs[14];// 14 -> N(12)
                    outputs[13] = temp2;     // 10 -> 0(13)
                    outputs[14] = inputs[15];// 15 -> P(14)
                    outputs[15] = temp;      // 9  -> Q(15)
                    if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 3ul)
                    {
                        Signal<T>::copy(ProcessorHarmonics<Hoa3d, T>::getNumberOfHarmonics() - 16ul, inputs+16, outputs+16);
                    }
                }
            }
        }

        //! This method number the channels from SID to ACN.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void numberFromSID(T const* inputs, T* outputs) noexcept
        {
            *(outputs++) = *(inputs++);
            for(size_t i = 1ul; i <= ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder(); i++)
            {
                for(size_t j = 0; j < 2ul * i + 1ul; j++)
                {
                    if(j%2)
                    {
                        m_harmonics[j] = inputs[size_t(j/2)];
                    }
                    else
                    {
                        m_harmonics[j] = inputs[2ul * i - size_t(j/2)];
                    }
                }
                Signal<T>::copy(i * 2ul + 1ul, m_harmonics, outputs);
                outputs += i * 2ul + 1ul;
                inputs += i * 2ul + 1ul;
            }
        }

        //! This method number the channels from ACN to SID.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void numberToSID(T const* inputs, T* outputs) noexcept
        {
            *(outputs++) = *(inputs++);
            for(size_t i = 1ul; i <= ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder(); i++)
            {
                for(size_t j = 0; j < 2ul * i + 1ul; j++)
                {
                    if(j%2)
                    {
                        m_harmonics[size_t(j/2)] = *(inputs++);
                    }
                    else
                    {
                        m_harmonics[2ul * i - size_t(j/2)] = *(inputs++);
                    }
                }
                Signal<T>::copy(i * 2ul + 1ul, m_harmonics, outputs);
                outputs += i * 2ul + 1ul;
            }
        }

        //! This method normalizes the channels from MaxN to SN3D.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void normalizeFromMaxN(T const* inputs, T* outputs) noexcept
        {
            outputs[0] = inputs[0] * T(sqrt(2.));
            outputs[1] = inputs[1] * T(sqrt(3.));
            outputs[2] = inputs[2] * T(sqrt(3.));
            outputs[3] = inputs[3] * T(sqrt(3.));
            if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 1ul)
            {
                outputs[4] = inputs[4] * T(sqrt(15.) / 2.);
                outputs[5] = inputs[5] * T(sqrt(15.) / 2.);
                outputs[6] = inputs[6] * T(sqrt(5.));
                outputs[7] = inputs[7] * T(sqrt(15.) / 2.);
                outputs[8] = inputs[8] * T(sqrt(15.) / 2.);
                if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 2ul)
                {
                    outputs[9]  = inputs[9]  * T(sqrt(35. / 8.));
                    outputs[10] = inputs[10] * T(sqrt(35.) / 3.);
                    outputs[11] = inputs[11] * T(sqrt(224. / 45));
                    outputs[12] = inputs[12] * T(sqrt(7.));
                    outputs[13] = inputs[13] * T(sqrt(224. / 45));
                    outputs[14] = inputs[14] * T(sqrt(35.) / 3.);
                    outputs[15] = inputs[15] * T(sqrt(35. / 8.));
                }
            }
            normalizeFromN3D(outputs, outputs);
        }

        //! This method normalizes the channels from SN3D to MaxN.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void normalizeToMaxN(T const* inputs, T* outputs) noexcept
        {
            outputs[0] = inputs[0] / T(sqrt(2.));
            outputs[1] = inputs[1] / T(sqrt(3.));
            outputs[2] = inputs[2] / T(sqrt(3.));
            outputs[3] = inputs[3] / T(sqrt(3.));
            if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 1ul)
            {
                outputs[4] = inputs[4] / T(sqrt(15.) / 2.);
                outputs[5] = inputs[5] / T(sqrt(15.) / 2.);
                outputs[6] = inputs[6] / T(sqrt(5.));
                outputs[7] = inputs[7] / T(sqrt(15.) / 2.);
                outputs[8] = inputs[8] / T(sqrt(15.) / 2.);
                if(ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder() > 2ul)
                {
                    outputs[9]  = inputs[9]  / T(sqrt(35. / 8.));
                    outputs[10] = inputs[10] / T(sqrt(35.) / 3.);
                    outputs[11] = inputs[11] / T(sqrt(224. / 45));
                    outputs[12] = inputs[12] / T(sqrt(7.));
                    outputs[13] = inputs[13] / T(sqrt(224. / 45));
                    outputs[14] = inputs[14] / T(sqrt(35.) / 3.);
                    outputs[15] = inputs[15] / T(sqrt(35. / 8.));
                }
            }
            normalizeToN3D(outputs, outputs);
        }

        //! This method normalizes the channels from N3D to SN3D.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void normalizeFromN3D(T const* inputs, T* outputs) noexcept
        {
            T norm = T(sqrt(3.));
            *(outputs++) = *(inputs++);
            *(outputs++) = *(inputs++) * norm;
            *(outputs++) = *(inputs++) * norm;
            *(outputs++) = *(inputs++) * norm;
            for(size_t i = 2; i <= ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder(); i++)
            {
                norm = T(sqrt(2. * T(i) + 1.));
                for(size_t j = 0; j < i * 2ul + 1ul; j++)
                {
                    *(outputs++) = *(inputs++) * norm;
                }
            }
        }

        //! This method normalizes the channels from SN3D to N3D.
        /**	You should use this method for in-place or not-in-place processing and sample by sample. The inputs array and outputs array contains the spherical harmonics samples and the minimum size must be the number of harmonics.
         @param     inputs   The inputs array.
         @param     outputs  The outputs array.
         */
        void normalizeToN3D(T const* inputs, T* outputs) noexcept
        {
            T norm = T(1. / sqrt(3.));
            *(outputs++) = *(inputs++);
            *(outputs++) = *(inputs++) * norm;
            *(outputs++) = *(inputs++) * norm;
            *(outputs++) = *(inputs++) * norm;
            for(size_t i = 2; i <= ProcessorHarmonics<Hoa3d, T>::getDecompositionOrder(); i++)
            {
                norm = T(1. / sqrt(2. * T(i) + 1.));
                for(size_t j = 0; j < i * 2ul + 1ul; j++)
                {
                    *(outputs++) = *(inputs++) * norm;
                }
            }
        }
        
        //! Retrieve the harmonic order of an input depending on the current numbering configuration.
        /** Retrieve the harmonic order of an input depending on the current numbering configuration.
         @param     index	The index of an harmonic.
         @return    The harmonic order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        long getInputHarmonicOrder(const size_t index) const noexcept
        {
            const Numbering numb = getNumbering();
            const bool acn = (numb == ACN) || (numb == toFurseMalham) || (numb == toSID);
            const bool malham = !acn && numb == fromFurseMalham;
            
            if(acn || (malham && index > 15))
            {
                return ProcessorHarmonics<Hoa3d, T>::getHarmonicOrder(index);
            }
            else
            {
                const long hdegree = long(sqrt(index));
                const long tmp = long(index) - (hdegree) * (hdegree);
                const long sign = (tmp % 2) ? -1l : 1l;
                return sign * long(hdegree + 1l - ((tmp - 1l) / 2. + 1l));
            }
        }
        
        //! Retrieve the harmonic order of an output depending on the current numbering configuration.
        /** Retrieve the harmonic order of an output depending on the current numbering configuration.
         @param     index	The index of an harmonic.
         @return    The harmonic order.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        long getOutputHarmonicOrder(const size_t index) const noexcept
        {
            const Numbering numb = getNumbering();
            const bool acn = (numb == ACN) || (numb == fromFurseMalham) || (numb == fromSID);
            const bool malham = !acn && numb == toFurseMalham;
            
            if(acn || (malham && index > 15))
            {
                return ProcessorHarmonics<Hoa3d, T>::getHarmonicOrder(index);
            }
            else
            {
                const long hdegree = long(sqrt(index));
                const long tmp = long(index) - (hdegree) * (hdegree);
                const long sign = (tmp % 2) ? -1l : 1l;
                return sign * long(hdegree + 1l - ((tmp - 1l) / 2. + 1l));
            }
        }
        
        //! Retrieve the name of an harmonic depending on the current numbering configuration.
        /** Retrieve the name of an harmonic depending on the current numbering configuration.
         This returns the name of the harmonic that contains its degree and its order for ACN and SID and harmonic name with letter code for FUMA.
         @param     index	The index of an harmonic.
         @param     isInput	Pass true to retrieve the input harmonic name, false for an output.
         @return    The name of the harmonic.
         @see       getHarmonicDegree()
         @see       getHarmonicOrder()
         */
        std::string getHarmonicName(const size_t index, const bool isInput) const noexcept
        {
            const Numbering numb = getNumbering();
            const bool acn = (numb == ACN) || (isInput && numb == toFurseMalham) || (isInput && numb == toSID) || (!isInput && numb == fromFurseMalham) || (!isInput && numb == fromSID);
            const bool malham = !acn && ((isInput && numb == fromFurseMalham) || (!isInput && numb == toFurseMalham));
            
            if(acn || (malham && index > 15))
            {
                // [0, 0], [1, -1], [1, 0], [1, 1], [2, -2], [2, -1], [2, 0], [2, 1], [2, 2] ...
                return ProcessorHarmonics<Hoa3d, T>::getHarmonicName(index);
            }
            else if(malham)
            {
                if(index < 16)
                {
                    static const char FurseMalham_3D_Letterings[16] =
                    {
                        'W', 'X', 'Y', 'Z',
                        'R', 'S', 'T', 'U', 'V',
                        'K', 'L', 'M', 'N', 'O', 'P', 'Q'
                    };
                    
                    std::string name = "Harmonic ";
                    return name += FurseMalham_3D_Letterings[index];
                }
                
                return "?";
            }
            else
            {
                // [0, 0], [1, 1], [1, -1], [1, 0], [2, 2], [2, -2], [2, 1], [2, -1], [2, 0] ...
                const long hdegree = long(sqrt(index));
                const long tmp = long(index) - (hdegree) * (hdegree);
                const long sign = (tmp % 2) ? -1l : 1l;
                const long horder = sign * long(hdegree + 1l - ((tmp - 1l) / 2. + 1l));
                std::ostringstream ostr;
                ostr <<  "Harmonic " << hdegree << " " << horder;
                return ostr.str();
            }
        }
    };

//! @endcond
}
