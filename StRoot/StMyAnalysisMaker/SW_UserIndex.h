#ifndef __SW_UserIndex__h__
#define __SW_UserIndex__h__

#include "FJ_includes.h"

FASTJET_BEGIN_NAMESPACE

class SW_UserIndex: public SelectorWorker {
    public:

    SW_UserIndex(int userindex): m_user_index(userindex) {};

    virtual bool pass(const PseudoJet & jet) const {
        // if (jet.user_index() == m_user_index)cout << Form("Asking to compare %i to %i ", jet.user_index(), m_user_index) << endl;
        return jet.user_index() == m_user_index;
    }

    /// returns a description of the worker
    virtual std::string description() const {
        std::ostringstream ostr;
        ostr << "User_Index Selected " << m_user_index << "\n";
        return ostr.str();
    }

    protected:
        int m_user_index;
};

inline Selector SelectorUserIndex (int i){
    return Selector(new SW_UserIndex(i));
}

FASTJET_END_NAMESPACE
#endif