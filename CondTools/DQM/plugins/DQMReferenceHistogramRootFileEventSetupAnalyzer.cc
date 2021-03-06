// C++ common header
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondFormats/Common/interface/FileBlob.h"
#include "CondFormats/DataRecord/interface/DQMReferenceHistogramRootFileRcd.h"
#include "DQMServices/Core/interface/DQMStore.h"

namespace edmtest {
  class DQMReferenceHistogramRootFileEventSetupAnalyzer : public edm::one::EDAnalyzer<> {
  public:
    typedef dqm::legacy::MonitorElement MonitorElement;
    typedef dqm::legacy::DQMStore DQMStore;
    explicit DQMReferenceHistogramRootFileEventSetupAnalyzer(const edm::ParameterSet& pset);
    explicit DQMReferenceHistogramRootFileEventSetupAnalyzer(int i);
    ~DQMReferenceHistogramRootFileEventSetupAnalyzer() override;
    void analyze(const edm::Event& event, const edm::EventSetup& setup) override;
    void beginRun(edm::Run const&, edm::EventSetup const&);

  private:
    const edm::ESGetToken<FileBlob, DQMReferenceHistogramRootFileRcd> fileBlobToken_;
    bool init_;
  };

  DQMReferenceHistogramRootFileEventSetupAnalyzer::DQMReferenceHistogramRootFileEventSetupAnalyzer(
      const edm::ParameterSet& ps)
      : fileBlobToken_(esConsumes<edm::Transition::BeginRun>()) {
    init_ = false;
    edm::LogPrint("DQMReferenceHistogramRootFileEventSetupAnalyzer")
        << "DQMReferenceHistogramRootFileEventSetupAnalyzer(const edm::ParameterSet &ps)" << std::endl;
  }

  DQMReferenceHistogramRootFileEventSetupAnalyzer::DQMReferenceHistogramRootFileEventSetupAnalyzer(int i)
      : fileBlobToken_(esConsumes<edm::Transition::BeginRun>()) {
    init_ = false;
    edm::LogPrint("DQMReferenceHistogramRootFileEventSetupAnalyzer")
        << "DQMReferenceHistogramRootFileEventSetupAnalyzer(int i) " << i << std::endl;
  }

  DQMReferenceHistogramRootFileEventSetupAnalyzer::~DQMReferenceHistogramRootFileEventSetupAnalyzer() {
    init_ = false;
    edm::LogPrint("DQMReferenceHistogramRootFileEventSetupAnalyzer")
        << "~DQMReferenceHistogramRootFileEventSetupAnalyzer" << std::endl;
  }

  void DQMReferenceHistogramRootFileEventSetupAnalyzer::analyze(const edm::Event& iEvent,
                                                                const edm::EventSetup& iSetup) {
    return;
  }

  void DQMReferenceHistogramRootFileEventSetupAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {
    edm::LogPrint("DQMReferenceHistogramRootFileEventSetupAnalyzer")
        << "DQMReferenceHistogramRootFileEventSetupAnalyzer::beginRun()" << std::endl;
    if (!init_) {
      init_ = true;
      edm::eventsetup::EventSetupRecordKey recordKey(
          edm::eventsetup::EventSetupRecordKey::TypeTag::findType("DQMReferenceHistogramRootFileRcd"));
      if (recordKey.type() == edm::eventsetup::EventSetupRecordKey::TypeTag()) {
        throw cms::Exception("Record not found") << "Record \"DQMReferenceHistogramRootFileRcd"
                                                 << "\" does not exist!" << std::endl;
      }
      const auto& rootgeo = &iSetup.getData(fileBlobToken_);
      edm::LogPrint("DQMReferenceHistogramRootFileEventSetupAnalyzer") << "ROOT FILE IN MEMORY" << std::endl;
      std::unique_ptr<std::vector<unsigned char> > tb((*rootgeo).getUncompressedBlob());
      // char filename[128];
      // sprintf(filename, "mem:%p,%ul", &(*tb)[0], (unsigned long) tb->size());
      // edm::Service<DQMStore>()->open(filename, false, "", "Reference");

      //here you can implement the stream for putting the TFile on disk...
      std::string outfile("dqmreference.root");
      std::ofstream output(outfile.c_str());
      output.write((const char*)&(*tb)[0], tb->size());
      output.close();

      DQMStore* dqm = &*edm::Service<DQMStore>();
      dqm->open(outfile, false, "", "Reference");
      remove(outfile.c_str());

      std::vector<MonitorElement*> mes = dqm->getAllContents("");
      for (std::vector<MonitorElement*>::iterator i = mes.begin(), e = mes.end(); i != e; ++i)
        edm::LogPrint("DQMReferenceHistogramRootFileEventSetupAnalyzer") << "ME '" << (*i)->getFullname() << "'\n";

      edm::LogPrint("DQMReferenceHistogramRootFileEventSetupAnalyzer") << "SIZE FILE = " << tb->size() << std::endl;
    }
  }

  DEFINE_FWK_MODULE(DQMReferenceHistogramRootFileEventSetupAnalyzer);
}  // namespace edmtest
