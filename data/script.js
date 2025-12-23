function loadContent(sectionId) {
  // 隐藏所有内容
  const sections = document.querySelectorAll('.section');
  sections.forEach(section => section.style.display = 'none');

  // 显示指定的内容
  const targetSection = document.getElementById(`${sectionId}-content`);
  if (targetSection) {
      targetSection.style.display = 'block';
  }
}

const contentMap = {
    'publications': 'publications-content',
    'xiangmu': 'xiangmu-content',
    'huojiang': 'huojiang-content',
    'about': 'about-content',
    'teach': 'teach-content',
  };

  // 显示全部内容
     function showAll() {
      for (let key in contentMap) {
        document.getElementById(contentMap[key]).style.display = 'block';
      }
    }

// 初始化时显示全部内容
document.addEventListener('DOMContentLoaded', function() {
  showAll();
});




