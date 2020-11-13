#ifndef LLPReco_DataFormats_FlatTableFiller_h
#define LLPReco_DataFormats_FlatTableFiller_h

#include "DataFormats/NanoAOD/interface/FlatTable.h"

template<class CLASS>
class FlatTableFiller
{
    public:
        virtual void push_back(const CLASS& obj) = 0;
        virtual void fill(std::unique_ptr<nanoaod::FlatTable>& flatTable) const = 0;
        
        virtual ~FlatTableFiller()
        {
        }
};

template<class CLASS>
class Property
{
    public:
        virtual const std::string& name() const = 0;
        virtual const std::string& doc() const = 0;
        virtual std::shared_ptr<FlatTableFiller<CLASS>> createFiller() const = 0;
        
        virtual ~Property()
        {
        }
};



template<class CLASS, class TYPE>
class PropertyTmpl:
    public Property<CLASS>
{
    private:
        TYPE CLASS::*memberPtr_;
        std::string name_;
        std::string doc_;
    public:
        PropertyTmpl(TYPE CLASS::*memberPtr, const std::string& name, const std::string& doc="doc"):
            memberPtr_(memberPtr),
            name_(name),
            doc_(doc)
        {
        }
        
        virtual ~PropertyTmpl()
        {
        }
        
        TYPE get(const CLASS& obj) const
        {
            return obj.*memberPtr_;
        }
    
        virtual const std::string& name() const
        {
            return name_;
        }
        
        virtual const std::string& doc() const
        {
            return doc_;
        }
        
        virtual std::shared_ptr<FlatTableFiller<CLASS>> createFiller() const;
};


template<class TYPE>
void fillTable(const std::vector<TYPE>& data, std::unique_ptr<nanoaod::FlatTable>& flatTable, const std::string& name, const std::string& doc)
{
    throw std::runtime_error(std::string("Filler not implemented for type '")+(typeid(TYPE)).name()+"'");
}

template<>
void fillTable<int>(const std::vector<int>& data, std::unique_ptr<nanoaod::FlatTable>& flatTable, const std::string& name, const std::string& doc)
{
    flatTable->addColumn<int>(name, data, doc, nanoaod::FlatTable::IntColumn);
}

template<>
void fillTable<float>(const std::vector<float>& data, std::unique_ptr<nanoaod::FlatTable>& flatTable, const std::string& name, const std::string& doc)
{
    flatTable->addColumn<float>(name, data, doc, nanoaod::FlatTable::FloatColumn);
}

template<class CLASS, class TYPE>
class FlatTableFillerTmpl:
    public FlatTableFiller<CLASS>
{
    private:
        PropertyTmpl<CLASS,TYPE> property_;
        std::vector<TYPE> data_;
    public:
        FlatTableFillerTmpl(const PropertyTmpl<CLASS,TYPE>& property):
            property_(property)
        {
        }
        
        virtual void push_back(const CLASS& obj)
        {
            data_.push_back(property_.get(obj));
        }
        
        virtual void fill(std::unique_ptr<nanoaod::FlatTable>& flatTable) const
        {
            fillTable<TYPE>(data_,flatTable,property_.name(),property_.doc());
        }

        virtual ~FlatTableFillerTmpl()
        {
        }
};



template<class CLASS, class TYPE> std::shared_ptr<FlatTableFiller<CLASS>> PropertyTmpl<CLASS,TYPE>::createFiller() const
{
    return std::shared_ptr<FlatTableFiller<CLASS>>(new FlatTableFillerTmpl<CLASS,TYPE>(*this));
}


template<class CLASS> using PropertyList = std::vector<std::shared_ptr<Property<CLASS>>>;

template<class CLASS>
class FlatTableFillerList
{
    protected:
        std::vector<std::shared_ptr<FlatTableFiller<CLASS>>> fillers_;
        
    public:
        FlatTableFillerList(const PropertyList<CLASS>& propertyList)
        {
            for (const auto& property: propertyList) fillers_.push_back(property->createFiller());
        }
        
        void push_back(const CLASS& obj)
        {
            for (const auto& filler: fillers_) filler->push_back(obj);
        }
        
        void fill(std::unique_ptr<nanoaod::FlatTable>& flatTable) const
        {
            for (const auto& filler: fillers_) filler->fill(flatTable);
        }
};

#define PROPERTY(class, name, doc) \
    std::shared_ptr<Property< class >>(new PropertyTmpl< class , decltype( class :: name )>(& class :: name , #name , doc))
   


#endif
